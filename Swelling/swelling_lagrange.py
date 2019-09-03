# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Based on the formualation in 2015 paper "Effect of solvent diffusion on
# crack-tip fields and driving force for fracture of hydrogels" which simplifies
# the free energy formulation in a prior 2015 paper, "A nonlinear, transient FE
# method for coupled solvent diffusion and large deformation of hydrogels"

# Import all packages in Dolfin module
from dolfin import *
from multiphenics import *

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 75,
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
        values[0] = (l0-1)*x[0]
        values[1] = (l0-1)*x[1]
        values[2] = (l0-1)*x[2]
        # Initial Chemical potential: mu0
        values[3] = 0 #ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
        # Initial lagrange multiplier, lambda
        values[4] = 0
    def value_shape(self):
         return (5,)

# Class for interfacing with the Newton solver
class Equation(NonlinearProblem):
    # Bilinear and linear forms are used to compute the residual and Jacobian
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)
        self.L = L              # Linear form
        self.a = a              # Bilinear form
        self.bcs = bcs          # Boundary conditions
    def F(self, b, x):          # Computes the residual vector "b"
        assemble(self.L, tensor=b)
        for bc in self.bcs:
            bc.apply(b)
    def J(self, A, x):          # Computes the Jacobian matrix "A"
        assemble(self.a, tensor=A)
        for bc in self.bcs:
            bc.apply(A)

# Class for marking everything on the domain except the top surface
class bulkWOsurf(SubDomain):
    # on_boundary will be false if you use pointwise method
    def inside(self, x, on_boundary):
        return x[1] <= 1.0 - tol

# Model parameters
#------------------------------------------------------------------------------
name = "penalty.xdmf"           # Name of file
tol = 1E-14                     # Tolerance for coordinate comparisons
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
l0 = 1.4                        # Initial Stretch (lambda_o)
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Global stepping and chemical stepping parameters
steps = 0                       # Steps (updated within loop)
c_steps = 0                     # Chemical step counter (updated within loop)
t_c_steps = 0                   # Total chemical steps
t_indent_steps = 5              # Total indentation steps
tot_steps = 5                   # Total number of time steps
# Indenter parameters
R = 0.25                        # Initial indenter radius
depth = 0.00                    # Initial indenter depth of indentation
# Time parameters
dt = 10**(-5)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression ("dt", dt=dt, degree=0)
t = 0.0                         # Initial time for paraview file
c_exp = 1.2                     # Controls time step increase (20%)

# Define mesh and mixed function space
#------------------------------------------------------------------------------
N_plane = 10                               # Number of elements on top plane
mesh = UnitCubeMesh(N_plane, 5, N_plane)   # Unit Cube
TT = TensorFunctionSpace(mesh,'DG',0)      # Tensor space for stress projection
V0 = FunctionSpace(mesh, "DG", 0)          # Vector space for contact pressure
# Taylor-Hood Elements for displacment (u) and chemical potential (mu)
#P2 = VectorElement("CG", mesh.ufl_cell(), 2)
#P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
P1 = FunctionSpace(mesh, "Lagrange", 1)
# Define mixed function space specifying underlying finite element
#V_element = BlockElement(P2, P1, P1)
#V = BlockFunctionSpace(mesh, V_element)
V = BlockFunctionSpace([P2, P1, P1], restrict=[None, None, bulkWOsurf()])

# Define functions in mixed function space V
#------------------------------------------------------------------------------
du = BlockTrialFunction(V)                       # Incremental trial function
v = BlockTestFunction(V)                         # Test Function
w = BlockFunction(V)                             # Current solution for u and mu
w0 = BlockFunction(V)                            # Previous solution for u and mu
# Split test functions and unknowns (produces a shallow copy not a deep copy)
(v_u, v_mu, v_l) = block_split(v)                 # Split test function
# Lambda is a reserved keyword in FEniCS
(u, mu, lmbda) = block_split(w)                   # Split current solution
(u0, mu0, lmbda0) = block_split(w0)               # Split previous solution

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

# Create mesh function over cell facets
#------------------------------------------------------------------------------
# An argument specifying the type of MeshFunction must be given (size_t). The
# second argument specifies the mesh, while the third argument gives the
# topological dimension (2 in this case)
facets = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
# Mark all facets as part of subdomain 0
facets.set_all(0)
# Top exterior facets are marked as subdomain 1, using 'top' boundary
top.mark(facets, 1)
# Measure redefines ds
ds = Measure('ds', subdomain_data=facets)

# Call class defining the whole domain except the top surface
bulkWOsurf = bulkWOsurf()
# Mark as part of subdomain 2
bulkWOsurf.mark(facets, 2)
# Visualization of marked facets
File("facets.pvd") << facets

# Bottom of the cube is attached to substrate and has fixed displacement from
# the reference relative to the current configuration
u_bot = Expression(("(l0-1)*x[0]", "0.0*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
    # Account for initial configuration: displacement in x & z direction not y

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(l0-1)*x[0]"), l0=l0, degree=1)
u_bf = Expression(("(l0-1)*x[2]"), l0=l0, degree=1)

# The Dirichlet BCs are specified in respective subspaces
# Displacement in first subspace V.sub(0)
bc_u = DirichletBC(V.sub(0), u_bot, bot)

# Roller displacement BCs on lateral faces
# Access normal degree of freedom (Ex. V.sub(0).sub(0) gives x direction)
bc_r_l = DirichletBC(V.sub(0).sub(0), u_lr, left)
bc_r_r = DirichletBC(V.sub(0).sub(0), u_lr, right)
bc_r_b = DirichletBC(V.sub(0).sub(2), u_bf, back)
bc_r_f = DirichletBC(V.sub(0).sub(2), u_bf, front)

# Test: Set lambda to 0 in domain
lmbda_val = Constant(0.0)
t1 = DirichletBC(V.sub(2), lmbda_val, bulkWOsurf, method="pointwise")

# Combined boundary conditions
bc = BlockDirichletBC([[bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f], [], [t1]])

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Initial conditions are created by using the class defined and then
# interpolating into a finite element space
init = InitialConditions(degree=1)          # Expression requires degree def.
w.BlockInterpolate(init)                        # Interpolate current solution
w0.interpolate(init)                       # Interpolate previous solution

# Kinematics
#------------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient from current time step
F0 = I + grad(u0)               # Deformation gradient from previous time step
CG = F.T*F                      # Right Cauchy-Green (CG) tensor

# Invariants of deformation tensors
Ic = tr(CG)                     # First invariant
J = det(F)                      # Current time step for third invariant
J0 = det(F0)                    # Previous time step for third invariant

# Definitions
#------------------------------------------------------------------------------
# Normalized nominal stress tensor: P = dU/dF
def P(u, mu):
    return F + (-1/J + (1/n)*(1/J + ln((J-1)/J) + chi/(J**2) - mu))*J*inv(F.T)

# Normalized flux
def Flux(u, mu):
    p1 = dot(inv(F), grad(mu))
    return -(J-1)*dot(p1,inv(F))

# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.

# Define the indenter with a parabola
indent = Expression("-d+(pow((x[0]),2)+pow((x[2]),2))/(2*R)", \
                        l0=l0, d=depth, R=R, degree=2)

# Variational problem
WF = [inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx \
    + lmbda*dot(v_u[1], ppos(u[1]-indent))*ds(1),
    (1/n)*((J-1)*v_mu*dx - (J0-1)*v_mu*dx - DT*dot(Flux(u, mu), grad(v_mu))*dx),
    v_l*ppos(u[1] - indent)*ds(1)]

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = block_derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = BlockNonlinearProblem(WF, w, bc, Jacobian)
solver_problem = BlockPETScSNESSolver(problem)
solver_problem.parameters.update(snes_solver_parameters["snes_solver"])

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile(name)

# Define contact pressure output (name for Paraview)
p = Function(V0, name="Contact pressure")

# Solve for each value using the previous solution as a starting point
while (steps < tot_steps):

    # Print to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))
    #print("Time: " + str(t))

    # Update fields containing u, mu, and lmbda and solve
    w0.block_vector()[:] = w.block_vector()
    solver_problem.solve()

    # Update total steps
    steps += 1
    # Update time parameters
    dt = dt*c_exp;                # Update time step with exponent value
    DT.dt = dt                    # Update time step for weak forms
    t += dt                       # Update total time for paraview file
    # Update indenter idndentation depth parameter
    if t_indent_steps > steps:
        depth += 0.01
        indent.depth = depth      # Update depth in expression class

    # This is a deep copy not a shallow copy like split(w) allowing us to write
    # the results to file
    (u, mu, lmbda) = w.block_split()
    PTensor = project(P(u, mu), TT)         # Project nominal stress

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential", "mu")
    PTensor.rename("Nominal Stress", "P")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u, t)
    file_results.write(mu, t)
    file_results.write(PTensor, t)
