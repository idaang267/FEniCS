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
# Class for marking everything on the domain except the top surface
class bulkWOsurf(SubDomain):
    # on_boundary will be false if you use pointwise method
    def inside(self, x, on_boundary):
        return x[1] <= 1.0 - tol

# Model parameters
#------------------------------------------------------------------------------
name = "Test_Lagrange_Hyperelasticity.xdmf"           # Name of file
tol = 1E-14                     # Tolerance for coordinate comparisons
# Elasticity parameters
E, nu = 10.0, 0.3                            # Young's Modulus, Poisson's ratio
mu = Constant(E/(2*(1 + nu)))                # Lame parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu))) # Lame parameter
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
start_point = Point(0, 0, 0)
end_point = Point(1.4, 1.4, 1.4)
# Cube
mesh = BoxMesh(start_point, end_point, N_plane, 4, N_plane)
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
V = BlockFunctionSpace([P2, P1], restrict=[None, bulkWOsurf()])

# Define functions in mixed function space V
#------------------------------------------------------------------------------
du = BlockTrialFunction(V)                       # Incremental trial function
v = BlockTestFunction(V)                         # Test Function
w = BlockFunction(V)                             # Current solution for u and mu
w0 = BlockFunction(V)                            # Previous solution for u and mu
# Split test functions and unknowns (produces a shallow copy not a deep copy)
(v_u, v_l) = block_split(v)                 # Split test function
# Lambda is a reserved keyword in FEniCS
(u, lmbda) = block_split(w)                   # Split current solution
(u0, lmbda0) = block_split(w0)               # Split previous solution

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube
# Note: these subdomains are set according to ParaView default view where
# x is right, y is up, and z-direction is out-of-plane
bot = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = l0)
# Lateral surfaces
left  = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = l0)
back  = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[2], side) && on_boundary", side = l0)

# Create mesh function over cell facets
#------------------------------------------------------------------------------
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
u_bot = Expression(("0.0*x[0]", "0.0*x[1]", "0.0*x[2]"), degree=1)
    # Account for initial configuration: displacement in x & z direction not y

# Roller boundary conditions (define normal directions) on lateral surfaces
u_roller = Expression(("0.0"), degree=1)

# The Dirichlet BCs are specified in respective subspaces
# Displacement in first subspace V.sub(0)
bc_u = DirichletBC(V.sub(0), u_bot, bot)

# Roller displacement BCs on lateral faces
# Access normal degree of freedom (Ex. V.sub(0).sub(0) gives x direction)
bc_r_l = DirichletBC(V.sub(0).sub(0), u_roller, left)
bc_r_r = DirichletBC(V.sub(0).sub(0), u_roller, right)
bc_r_b = DirichletBC(V.sub(0).sub(2), u_roller, back)
bc_r_f = DirichletBC(V.sub(0).sub(2), u_roller, front)

# Test: Set lambda to 0 in domain
lmbda_val = Constant(0.0)
t1 = DirichletBC(V.sub(1), lmbda_val, bulkWOsurf, method="pointwise")

# Combined boundary conditions
bc = BlockDirichletBC([[bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f], [t1]])

# Kinematics
#------------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient from current time step
CG = F.T*F                      # Right Cauchy-Green (CG) tensor

# Invariants of deformation tensors
Ic = tr(CG)                     # First invariant
J = det(F)                      # Current time step for third invariant

# Definitions
#------------------------------------------------------------------------------
# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.

# Define the indenter with a parabola
indent = Expression("-d+(pow((x[0]),2)+pow((x[2]),2))/(2*R)", \
                        l0=l0, d=depth, R=R, degree=2)

# Define the strain energy density and the total potential energy
#------------------------------------------------------------------------------
# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3 - 2*ln(J)) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = [psi*dx,
      lmbda*dot(v_u[1], ppos(u[1]-indent))*ds(1) + v_l*ppos(u[1] - indent)*ds(1)]

# Compute first variation of Pi
WF = block_derivative(Pi, w, v)
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
    (u, lmbda) = w.block_split()

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    PTensor.rename("Nominal Stress", "P")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u, t)
