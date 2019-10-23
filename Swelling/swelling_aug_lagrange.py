# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Based on the formualation in 2015 paper "Effect of solvent diffusion on
# crack-tip fields and driving force for fracture of hydrogels" which simplifies
# the free energy formulation in a prior 2015 paper, "A nonlinear, transient FE
# method for coupled solvent diffusion and large deformation of hydrogels"

# Import modules
from dolfin import *        # Dolfin Module
from multiphenics import *  # For restricting a variable to a certain domain
from mshr import *          # For generating domains

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 200,
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

# Model parameters
#------------------------------------------------------------------------------
name = "AugL.xdmf"              # Name of file
k_pen = 1E4                     # Penalty Parameter
l0 = 1.4                        # Initial Stretch (lambda_o)
# Mesh Resolution
inPlaneRes = 16                 # Along top surface (without refinement)
outPlaneRes = 8                 # Resolution along cube height
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Stepping parameters
steps = 0                       # Steps (updated within loop)
tot_steps = 10                  # Total number of time steps for indentation
# Indenter parameters
R = 0.25                        # Initial indenter radius
depth = 0.00                    # Initial indenter depth of indentation
depth_indent = 0.01             # Amount to indent in each iteration
# Time parameters
dt = 10**(-5)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression ("dt", dt=dt, degree=0)
t = 0.0                         # Initial time for paraview file
c_exp = 1.1                     # Controls time step increase (20%)

# Create specific boundaries using classes
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube. These subdomains are
# set according to ParaView default view where x is right, y is up, and
# z-direction is out-of-plane

# Bottom and top
bot = CompiledSubDomain("near(x[1], 0.0) && on_boundary")
class top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary
# Lateral surfaces
class symmetry_x(SubDomain):    # left and right surfaces
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) or near(x[0], 1.0) and on_boundary
class symmetry_z(SubDomain):    # front and back surfaces
    def inside(self, x, on_boundary):
        return near(x[2], 0.0) or near(x[2], 1.0) and on_boundary
# Region where mesh is refined
class refDomain(SubDomain):
    def inside(self, x, on_boundary):
        return (pow(x[0]-0.5,2) <= pow(0.35,2) - pow(x[1]-1,2) - pow(x[2]-0.5,2) - DOLFIN_EPS)

# Convert classes
sym_x = symmetry_x()
sym_z = symmetry_z()
top = top()
refDomain = refDomain()

# Create and define mesh using mshr
#------------------------------------------------------------------------------
# Create unit cube mesh
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), inPlaneRes, outPlaneRes, inPlaneRes)

# Refinement of domain using the refine() function
d = mesh.topology().dim()           # Topology dimension
r_markers = MeshFunction("bool", mesh, d, False)
refDomain.mark(r_markers, True)
refinedMesh = refine(mesh,r_markers)

# Transfering a non-negative integer-valued (`size_t`) `MeshFunction` to the
# refined mesh using the `adapt()` function:
meshFunctionToAdapt = MeshFunction("size_t", mesh, d, 0)
refDomain.mark(meshFunctionToAdapt,1)
adaptedMeshFunction = adapt(meshFunctionToAdapt,refinedMesh)

# Create subdomains from refined mesh
subdomains = MeshFunction("size_t", refinedMesh, refinedMesh.topology().dim(), refinedMesh.domains())
# Create boundaries from refined mesh
boundaries = MeshFunction("size_t", refinedMesh, refinedMesh.topology().dim() - 1, 0)

# Marking for visualization
boundaries.set_all(0)
top.mark(boundaries, 1)

# Set interface restriction
top_restrict = MeshRestriction(refinedMesh, top)

# Measure redefines dx and ds
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Define mixed function space
#------------------------------------------------------------------------------
# Function spaces for projection
TT = TensorFunctionSpace(refinedMesh,'DG',0)    # Tensor space for stress
V0 = FunctionSpace(refinedMesh, "DG", 0)        # Vector space for contact pressure
V2 = FunctionSpace(refinedMesh, "Lagrange", 1)  # For Gap function
# Taylor-Hood Elements for displacment (u), chemical potential (mu), lambda (lmbda)
P2 = VectorFunctionSpace(refinedMesh, "Lagrange", 2)
P1 = FunctionSpace(refinedMesh, "Lagrange", 1)
# Define mixed function space where lambda is restricted to the top surface
# and should remain 0 for the bulk
V = BlockFunctionSpace([P2, P1, P1], restrict=[None, None, top_restrict])

# Define functions in mixed function space V
#------------------------------------------------------------------------------
du = BlockTrialFunction(V)                  # Incremental trial function
v = BlockTestFunction(V)                    # Test Function
w = BlockFunction(V)                        # Current solution for u and mu
w0 = BlockFunction(V)                       # Previous solution for u and mu
# Split test functions and unknowns (produces a shallow copy not a deep copy)
(v_u, v_mu, v_l) = block_split(v)           # Split test function
# Lambda is a reserved keyword in FEniCS
(u, mu, lmbda) = block_split(w)             # Split current solution
(u0, mu0, lmbda0) = block_split(w0)         # Split previous solution

# Define gap and augmented lmbda output (name for Paraview)
gap = Function(V2, name="Gap")
aug_lmbda = Function(P1, name="Augmented Lagrangian")

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Bottom of the cube is attached to substrate and has fixed displacement from
# the reference relative to the current configuration
u_bot = Expression(("(l0-1)*x[0]", "0.0*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
    # Account for initial configuration: displacement in x & z direction not y

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(l0-1)*x[0]"), l0=l0, degree=1)
u_bf = Expression(("(l0-1)*x[2]"), l0=l0, degree=1)

# The Dirichlet BCs are specified in respective subspaces
bc_u = DirichletBC(V.sub(0), u_bot, bot)    # u in first subspace V.sub(0)

# Roller displacement BCs on lateral faces
# Access normal degree of freedom (Ex. V.sub(0).sub(0) gives x direction)
bc_r_l = DirichletBC(V.sub(0).sub(0), u_lr, sym_x)
bc_r_r = DirichletBC(V.sub(0).sub(0), u_lr, sym_x)
bc_r_b = DirichletBC(V.sub(0).sub(2), u_bf, sym_z)
bc_r_f = DirichletBC(V.sub(0).sub(2), u_bf, sym_z)

# Combined boundary conditions for each subspace
bc = BlockDirichletBC([[bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f], [], []])

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Extract subfunction
block_u = w.sub(0)
block_p = w.sub(1)
block_l = w.sub(2)
# Interpolate on subfunction
u_ini = Expression(("(l0-1)*x[0]", "(l0-1)*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
block_u.interpolate(u_ini)
block_p.interpolate(Constant(0.0))
block_l.interpolate(Constant(0.0))
# Update the original block functions by assign(receiving_funcs, assigning_func)
# Current solution
assign(w.sub(0), block_u)
assign(w.sub(1), block_p)
assign(w.sub(2), block_l)
# Previous solution
assign(w0.sub(0), block_u)
assign(w0.sub(1), block_p)
assign(w0.sub(2), block_l)

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
# Define the augmented lagrangian
def aug_l(x):
    return lmbda + k_pen*(indent - u[1])

# Define the indenter with a pabola
indent = Expression("-depth+(l0-1)+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)", \
                        l0=l0, depth=depth, R=R, degree=2)

# Variational problem
WF = [inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx \
        - aug_l(lmbda)*v_u[1]*ds + ppos(aug_l(lmbda))*v_u[1]*ds,
    (1/n)*((J-1)*v_mu*dx - (J0-1)*v_mu*dx - DT*dot(Flux(u, mu), grad(v_mu))*dx),
    (indent-u[1])*v_l*ds - (1/k_pen)*ppos(aug_l(lmbda))*v_l*ds]

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = block_derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = BlockNonlinearProblem(WF, w, bc, Jacobian)
solver = BlockPETScSNESSolver(problem)
solver.parameters.update(snes_solver_parameters["snes_solver"])
# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile(name)

# Solve for each value using the previous solution as a starting point
while (steps <= tot_steps):

    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))
    print("Time: " + str(t))

    # Update fields containing u, mu, and lmbda and solve
    w0.block_vector()[:] = w.block_vector()
    solver.solve()

    # Update total steps
    steps += 1
    # Update time parameters
    dt = dt*c_exp;                  # Update time step with exponent value
    DT.dt = dt                      # Update time step for weak forms
    t += dt                         # Update total time for paraview file
    # Update indenter depth parameter then update in expression class
    depth += depth_indent
    indent.depth = depth

    # Write the results to a file using a deep copy not a shallow copy
    (u, mu, lmbda) = w.block_split()
    PTensor = project(P(u, mu), TT)                 # Project nominal stress
    gap.assign(project(indent - u[1], V2))          # Project gap into V2 space
    aug_lmbda.assign(project(aug_l(lmbda), P1))     # Project augmented lambda
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential", "mu")
    gap.assign(project(indent - u[1], V2))        # Project gap into V2 space
    lmbda.rename("Lagrange Multiplier", "lmbda")
    PTensor.rename("Nominal Stress", "P")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u, t)
    file_results.write(mu, t)
    file_results.write(lmbda, t)
    file_results.write(gap, t)
    file_results.write(PTensor, t)
    file_results.write(aug_lmbda, t)
