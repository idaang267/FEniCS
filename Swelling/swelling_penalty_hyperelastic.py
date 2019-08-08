# First, the required modules are imported
from dolfin import *
import matplotlib.pyplot as plt

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 100,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          # Newton line search
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "error_on_nonconvergence": False}}

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4

# Defining Classes
#------------------------------------------------------------------------------
# Initial condition (IC) class for displacement and chemical potential
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1], values[2])
        values[0] = (0.0)*x[0]
        values[1] = (0.0)*x[1]
        values[2] = (0.0)*x[2]
    def value_shape(self):
         return (3,)

# Define all parameters
#------------------------------------------------------------------------------
B  = Constant((0.0, 0.0, 0.0))               # Body force per unit volume
# Elasticity parameters
E, nu = 10.0, 0.3                            # Young's Modulus, Poisson's ratio
mu = Constant(E/(2*(1 + nu)))                # Lame parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu))) # Lame parameter
#
steps = 0                       # Steps (updated within loop)
t_steps = 50                     # Total number of time steps
# Indenter parameters
R = 0.25                        # Indenter radius
depth = 0.01                   # Indenter depth
l0 = 1.4                        # Initial Stretch (lambda_o)
# Penalty parameter
pen = Constant(1e4)
# Time parameters
dt = 10**(-5)                   # Starting time step
# NOTE: Time step must be in expression for updating the weak form
t = 0.0                        # Initial time for paraview file
c_exp = 1.2                    # Controls time step increase (20%)

# Mesh and function spaces
#------------------------------------------------------------------------------
# Create tetrahedral mesh of the domain and a function space on this mesh
mesh = UnitCubeMesh(8, 3, 8)
# Vector space for contact pressure
V0 = FunctionSpace(mesh, "DG", 0)
# Function space for displacement (u)
V = VectorFunctionSpace(mesh, "Lagrange", 2)

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube
# Note: these subdomains are set according to ParaView default view where
# x is right, y is up, and z-direction is out-of-plane
bot = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
# Lateral surfaces
left  = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
back  = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)

# Exterior facets
#------------------------------------------------------------------------------
facets = MeshFunction("size_t", mesh, 2)
# Mark all facets as part of subdomain 0
facets.set_all(0)
# Top exterior facets are marked as subdomain 1, using 'top' boundary
top.mark(facets, 1)
# Measure redefines ds
ds = Measure('ds', subdomain_data=facets)

# Bottom of the cube is attached to substrate and has fixed displacement
u_bot = Expression(("0.0*x[0]", "0.0*x[1]", "0.0*x[2]"), degree=1)

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("0.0*x[0]"), degree=1)
u_bf = Expression(("0.0*x[2]"), degree=1)

# Displacement defined for function space V
bc_u = DirichletBC(V, u_bot, bot)

# Roller displacement BCs on lateral faces, access normal degree of freedom
# V.sub(0).sub(0) gives x direction)
bc_r_l = DirichletBC(V.sub(0), u_lr, left)
bc_r_r = DirichletBC(V.sub(0), u_lr, right)
# V.sub(0).sub(2) gives z direction
bc_r_b = DirichletBC(V.sub(2), u_bf, back)
bc_r_f = DirichletBC(V.sub(2), u_bf, front)

# Combined boundary conditions
bcs = [bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f]

# Trial and test functions, and the most recent approximate displacement, u are
# defined on the finite element space V.
# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
u0 = Function(V)

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Initial conditions are created by using the class defined and then
# interpolating into a finite element space
init = InitialConditions(degree=1)          # Expression requires degree def.
u.interpolate(init)                         # Interpolate current solution
u0.interpolate(init)                        # Interpolate previous solution

# With the functions defined, the kinematic quantities involved in the model
# are defined using UFL syntax

# Kinematics
d = len(u)          # Length of displacement vector
I = Identity(d)     # Identity tensor
F = I + grad(u)     # Deformation gradient
C = F.T*F           # Right Cauchy-Green tensor C = F^T F

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Set material parameters and define the strain energy density and the total
# potential energy

# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.

# Indenter parameters
# Define the indenter with a parabola
indent = Expression("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)", \
                        depth=depth, R=R, degree=2)

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - dot(B, u)*dx + pen*dot(u[1], ppos(u[1]-indent))*ds(1)

# Directional derivatives are now computed of Pi and L

# Compute first variation of Pi
F = derivative(Pi, u, v) # directional derivative about u in the direction of v
# Compute Jacobian of F
Jacobian = derivative(F, u, du)

# The complete variational problem can now be solved by a single call to solve
problem = NonlinearVariationalProblem(F, u, bcs, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

# Solution 'u' is saved to a file in xdmf format,
file_results = XDMFFile("hyperelastic_penalty.xdmf"); # Save solution in VTK format

# Define contact pressure output (name for Paraview)
p = Function(V0, name="Contact pressure")

while (steps < t_steps):
    # Print to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))
    print("Time: " + str(t))
    # Solve for updated indentation
    u0.vector()[:] = u.vector()
    solver_problem.solve()

    # Update total steps
    steps += 1
    # Update indenter indentation depth parameter
    depth += 0.01
    indent.depth = depth          # Update depth in expression class
    # Update time parameters
    dt = dt*c_exp;                # Update time step with exponent value
    t += dt                       # Update total time for paraview file

    p.assign(project(pen*ppos(u[1]-indent), V0))    # Project pressure

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u, t)
    file_results.write(p, t)
