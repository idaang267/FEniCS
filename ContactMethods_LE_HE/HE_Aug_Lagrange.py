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

# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          'absolute_tolerance':1e-5,
                                          'relative_tolerance':1e-5,
                                          "maximum_iterations": 200,
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

# Parameters
# -----------------------------------------------------------------------------
name = "HE_AugL.xdmf"       # File name
E = 10.0                    # Material Parameter: Young's Modulus
nu = 0.3                    # Material Parameter: Poisson's Ratio
# Lame parameters
mu = Constant(E/(2*(1 + nu)))
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))
# Definitions of body force and traction
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
# Stepping parameters
steps = 0
tot_steps = 40
# depth
depth = 0.00
R = 0.25
# Augmented lagrangian constant: epsilon_N
k_pen = 1E4

# Create Mesh and define subdomains for restrictions and boundary conditions
# -----------------------------------------------------------------------------
domain = Box(Point(0., 0., 0.), Point(1.0, 1.0, 1.0))
mesh = generate_mesh(domain, 15)

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
V2 = FunctionSpace(mesh, "CG", 1)       # To plot Gap function
# Block mixed equal order function space
W = BlockFunctionSpace([V, V1], restrict=[None, boundary_restriction])

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, lagrange multiplier)
du = BlockTrialFunction(W)
ul = BlockFunction(W)
v = BlockTestFunction(W)
gap = Function(V2, name="Gap")

# Split trial and test functions into displacement and lagrange multiplier
(u, l) = block_split(ul)
(v_u, v_l) = block_split(v)

# Boundary Conditions
# -----------------------------------------------------------------------------
# Define Dirichlet boundary (x = 0 or x = 1)
c = Constant((0.0, 0.0, 0.0))
# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(0.0)*x[0]"), degree=1)
u_bf = Expression(("(0.0)*x[2]"), degree=1)

# Lateral surfaces
class symmetry_x(SubDomain):    # left and right surfaces
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) or near(x[0], 1.0) and on_boundary
class symmetry_z(SubDomain):    # front and back surfaces
    def inside(self, x, on_boundary):
        return near(x[2], 0.0) or near(x[2], 1.0) and on_boundary

sym_x = symmetry_x()
sym_z = symmetry_z()

# The Dirichlet BCs are specified in a subspace
bc_bot = DirichletBC(W.sub(0), c, bot)
# Sub displacement subspace to obtain direction of interest
bc_r_l = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_r = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_b = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)
bc_r_f = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)

bcs = BlockDirichletBC([[bc_bot, bc_r_l, bc_r_r, bc_r_b, bc_r_f], []])

# Kinematics
# -----------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor
J  = det(F)                     # 3rd invariant of the deformation tensor

# Constitutive relation
# -----------------------------------------------------------------------------
# P = dW/dF:
def P(u):
    return mu*(F - inv(F.T)) + lmbda*ln(J)*inv(F.T)
# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.
# Define the augmented lagrangian
def aug_l(x):
    return l + k_pen*(u[1] - indent)

# Expression for indenter
indent = Expression(("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)"), \
                    depth=depth, R=R, degree=0)

# Weak Form
F = [inner(P(u), grad(v_u))*dx - dot(B, v_u)*dx - dot(T, v_u)*ds + v_u[1]*ppos(aug_l(l))*ds,
    v_l*(u[1]-indent)*ds - (1/k_pen)*l*v_l*ds]

# Directional derivative
Jacobian = block_derivative(F, ul, du)

# Solve variational problem
varproblem = BlockNonlinearProblem(F, ul, bcs, Jacobian)
solver = BlockPETScSNESSolver(varproblem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

# Save results to an .xdmf file since we have multiple fields
file_results = XDMFFile(name)

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
while steps <= tot_steps:
    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))

    # Solve the nonlinear problem
    solver.solve()

    # Update the displacement value
    depth += 0.005
    indent.depth = depth
    steps += 1
    # Split results
    (u, l) = ul.block_split()
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    l.rename("Lagrange Multiplier", "l")
    gap.assign(project(indent-u[1], V2))
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u,steps)
    file_results.write(l, steps)
    file_results.write(gap, steps)
