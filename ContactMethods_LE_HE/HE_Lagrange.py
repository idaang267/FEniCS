# Hyperelasticity

# Using the Compressible Hyperelasticity Strain Energy Formulation
#   The stored strain energy density:
#   Compressible neo-Hookean model:
#       W = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2
#   Incompressible neo-Hookean model (J = det(F) = 1):
#       W = (mu/2)*(Ic-3)

# Edited to use SNES solver

# Import modules
from dolfin import *
from multiphenics import *
from mshr import *
#import numpy as np

# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          'absolute_tolerance':1e-5,
                                          'relative_tolerance':1e-5,
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"
# The following two fields were in the original hyperelasticity demo
parameters["form_compiler"]["cpp_optimize"] = True
# For processing UFLACS - unified form language expressions
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["allow_extrapolation"] = True

solver_par = NonlinearVariationalSolver.default_parameters()
solver_par.rename("solver")

# Material Parameters
E = 10.0
nu = 0.3
# Elasticity parameters
mu = Constant(E/(2*(1 + nu)))
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))
# Definitions of body force and traction
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
# Stepping parameters
steps = 0
tot_steps = 6
# depth
depth = 0.00
R = 0.25
# Augmented lagrangian constant: epsilon_N
k_pen = 1E4

# Create Mesh
domain = Box(Point(0., 0., 0.), Point(1.0, 1.0, 1.0))
mesh = generate_mesh(domain, 10)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
# Create boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Define class
class top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary

class bot(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

# Mark boundary to restrict the lagrange multiplier to
top = top()
bot = bot()
top.mark(boundaries, 1)

# Dirichlet boundary
boundary_restriction = MeshRestriction(mesh, top)

# Measure
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Taylor-Hood Function space
V = VectorFunctionSpace(mesh, "Lagrange", 2)
V1 = FunctionSpace(mesh, "Lagrange", 2)
# Block mixed equal order function space
W = BlockFunctionSpace([V, V1], restrict=[None, boundary_restriction])

# Define Dirichlet boundary (x = 0 or x = 1)
c = Constant((0.0, 0.0, 0.0))
# The Dirichlet BCs are specified in a subspace
bc_bot = DirichletBC(W.sub(0), c, bot)
bcs = BlockDirichletBC([[bc_bot],[]])

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, lagrange multiplier)
# Define trial and test functions
du = BlockTrialFunction(W)
ul = BlockFunction(W)
v = BlockTestFunction(W)
# Split trial and test functions into displacement and lagrange multiplier
(u, l) = block_split(ul)
(v_u, v_l) = block_split(v)

# Kinematics
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor
J  = det(F)                     # 3rd invariant of the deformation tensor

# Where P = dW/dF:
def P(u):
    return mu*(F - inv(F.T)) + lmbda*ln(J)*inv(F.T)
# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.
# Expression for indenter
indent = Expression(("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)"), \
                    depth=depth, R=R, degree=0)

# Expression for indenter
indent = Expression(("-depth"), \
                    depth=depth, R=R, degree=0)

# Weak Form
F = [inner(P(u), grad(v_u))*dx - dot(B, v_u)*dx - dot(T, v_u)*ds + inner(l,v_u[1])*ds + inner(l+k_pen*(u[1]-indent),v_u[1])*ds,
    v_l*(u[1]-indent)*ds - (1/k_pen)*l*v_l*ds]

# Directional derivative
Jacobian = block_derivative(F, ul, du)

# Solve variational problem
varproblem = BlockNonlinearProblem(F, ul, bcs, Jacobian)
solver = BlockPETScSNESSolver(varproblem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

# Save results to an .xdmf file since we have multiple fields
file_results = XDMFFile("Lagrange_HE.xdmf")

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
while steps <= tot_steps:

    # Solve the nonlinear problem (using Newton solver)
    solver.solve()

    # Update the displacement value
    depth += 0.01
    indent.depth = depth
    steps += 1
    # Split results
    (u, l) = ul.block_split()
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    l.rename("Lagrange Multiplier", "l")

    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u,steps)
    file_results.write(l, steps)
