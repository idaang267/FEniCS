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
                                          'absolute_tolerance': 1e-5,
                                          'relative_tolerance': 1e-5,
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
name = "Test.xdmf"               # File name
# Mesh resolution parameters
inPlaneRes = 15
outPlaneRes = 8
# Material Parameters
E = 10.0                            # Young's Modulus
nu = 0.3                            # Poisson's Ratio
mu = E/(2*(1 + nu))                 # First Lame parameters
lmbda = E*nu/((1 + nu)*(1 - 2*nu))  # Second Lame parameter
# Definitions of body force and traction
B  = Constant((0.0, 0.0, 0.0))      # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))      # Traction force on the boundary
# Stepping parameters
steps = 0                           # Updates in the loop
tot_steps = 20                      # Total amount of steps
# Indenter parameters
depth = 0.00                        # Initial Depth (Updates in loop)
depth_update = 0.01                 # Update of initial depth in each iteration
R = 0.25                            # Fixed radius of indenter
# Augmented lagrangian constant
k_pen = 1E3

# Define classes
# -----------------------------------------------------------------------------
class top(SubDomain):           # Top surface
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary
class bot(SubDomain):           # Bottom surface
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

# Lateral surfaces
class symmetry_x(SubDomain):    # left and right surfaces
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) or near(x[0], 1.0) and on_boundary
class symmetry_z(SubDomain):    # front and back surfaces
    def inside(self, x, on_boundary):
        return near(x[2], 0.0) or near(x[2], 1.0) and on_boundary

class refDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (pow(x[0]-0.5,2) <= pow(0.25,2) - pow(x[1]-1,2) - pow(x[2]-0.5,2) - DOLFIN_EPS)

# Convert classes
top = top()
bot = bot()
sym_x = symmetry_x()
sym_z = symmetry_z()
refDomain = refDomain()

# Create Mesh and define subdomains for restrictions and boundary conditions
# -----------------------------------------------------------------------------
# Define mesh and subdomain:
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), inPlaneRes, outPlaneRes, inPlaneRes)
d = mesh.topology().dim()

# Refinement using the `refine()` function and a Boolean `MeshFunction`:
r_markers = MeshFunction("bool", mesh, d, False)
refDomain.mark(r_markers, True)
refinedMesh = refine(mesh, r_markers)

# Transfering a non-negative integer-valued (`size_t`) `MeshFunction` to the
# refined mesh using the `adapt()` function:
meshFunctionToAdapt = MeshFunction("size_t", mesh, d, 0)
refDomain.mark(meshFunctionToAdapt,1)
adaptedMeshFunction = adapt(meshFunctionToAdapt,refinedMesh)

# Create subdomains
subdomains = MeshFunction("size_t", refinedMesh, refinedMesh.topology().dim(), refinedMesh.domains())
# Create boundaries
boundaries = MeshFunction("size_t", refinedMesh, refinedMesh.topology().dim() - 1)

# Mark top boundary to restrict the lagrange multiplier
top.mark(boundaries, 1)
boundary_restriction = MeshRestriction(refinedMesh, top)

# Measure redefines the subdomains and boundaries
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Define all function spaces
# -----------------------------------------------------------------------------
# Define Taylor-Hood function spaces
V = VectorFunctionSpace(refinedMesh, "Lagrange", 2)
V1 = FunctionSpace(refinedMesh, "Lagrange", 2)
# Block mixed equal order function space
W = BlockFunctionSpace([V, V1], restrict=[None, boundary_restriction])
# Function spaces for interpolation
V2 = FunctionSpace(refinedMesh, "Lagrange", 1)       # For Gap function

# Define trial and test functions and unknown solutions
du = BlockTrialFunction(W)
ul = BlockFunction(W)           # (u,l) = (displacement, lagrange multiplier)
v = BlockTestFunction(W)
# Create gap function and augmented lagrangian (name for paraview)
gap = Function(V2, name="Gap")
lagrange = Function(V2, name="Augmented Lagrangian")

# Split trial and test functions into displacement and lagrange multiplier
(u, l) = block_split(ul)
(v_u, v_l) = block_split(v)

# Boundary Conditions
# -----------------------------------------------------------------------------
# Define Dirichlet boundary
u_bot = Expression(("0.0*x[0]", "0.0*x[1]", "0.0*x[2]"),  degree=1)
# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("0.0*x[0]"), degree=1)
u_bf = Expression(("0.0*x[2]"), degree=1)

# The Dirichlet BCs are specified in the first subspace
# Fixed bottom boundary condition
bc_bot = DirichletBC(W.sub(0), u_bot, bot)
# Sub displacement subspace to obtain direction of interest
bc_r_l = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_r = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_b = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)
bc_r_f = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)

# Combine all boundary conditions (in this case all displacement bcs)
bcs = BlockDirichletBC([[bc_bot, bc_r_l, bc_r_r, bc_r_b, bc_r_f], []])

# Kinematics
# -----------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor
J = det(F)                      # 3rd invariant of the deformation tensor

# Constitutive relation
# -----------------------------------------------------------------------------
def P(u):                       # P = dW/dF:
    return mu*(F - inv(F.T)) + lmbda*ln(J)*inv(F.T)
def ppos(x):                    # Define the Mackauley bracket <x>+
    return (x+abs(x))/2.
def aug_l(x):                   # Define the augmented lagrangian
    return l + k_pen*(indent-u[1])

# Parabolic expression for indenter
indent = Expression(("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)"), \
                    depth=depth, R=R, degree=0)

# Two equations for the weak form
F = [inner(P(u), grad(v_u))*dx - dot(B, v_u)*dx - dot(T, v_u)*ds \
    - aug_l(l)*v_u[1]*ds + ppos(aug_l(l))*v_u[1]*ds,
    (indent-u[1])*v_l*ds - (1/k_pen)*ppos(aug_l(l))*v_l*ds]

# Directional derivative
Jacobian = block_derivative(F, ul, du)

# Solve variational problem using PETSc SNES solver
varproblem = BlockNonlinearProblem(F, ul, bcs, Jacobian)
solver = BlockPETScSNESSolver(varproblem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

# Save results to an .xdmf file since we have multiple fields
file_results = XDMFFile(name)

# Solve with Newton solver using the previous solution as a starting point
while steps <= tot_steps:
    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))
    # Solve the nonlinear problem
    solver.solve()
    steps += 1                  # Update steps
    depth += depth_update       # Update the depth of the indenter
    indent.depth = depth        # Update depth in expression class
    (u, l) = ul.block_split()   # Split results (hard not soft copy)
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    l.rename("Lagrange Multiplier", "l")
    gap.assign(project(indent - u[1], V2))        # Project gap into V2 space
    lagrange.assign(project(aug_l(l), V2))
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u, steps)
    file_results.write(l, steps)
    file_results.write(gap, steps)
    file_results.write(lagrange, steps)
