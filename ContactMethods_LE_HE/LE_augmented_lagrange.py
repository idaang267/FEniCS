# Increasing complexity

# Import all modules
from dolfin import *
from multiphenics import *
from mshr import *
from ufl import rank

snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": False}}

# Parameters
name = "LE_Aug_Lagrange.xdmf"
# Steps
steps = 0
tot_steps = 20
# Define the indenter (obstacle) with a parabola
R = 0.25                            # Radius
depth = 0.00                        # Depth
# Material Parameters
E = Constant(10.0)                  # Young's Modulus
nu = Constant(0.3)                  # Poisson's Ratio
mu = E/(2*(1+nu))                   # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))      # Lamé Coefficient
# Initial Stretch (lambda_o)
l0 = 1.4
# Augmented lagrangian constant: epsilon_N
k_pen = 1E4

# Create Mesh
domain = Box(Point(0., 0., 0.), Point(1.0, 1.0, 1.0))
mesh = generate_mesh(domain, 25)

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

# Function space
V = VectorFunctionSpace(mesh, "Lagrange", 2)
V1 = FunctionSpace(mesh, "Lagrange", 2)
V2 = FunctionSpace(mesh, "CG", 1) # Gap function
# Block mixed equal order function space
W = BlockFunctionSpace([V, V1], restrict=[None, boundary_restriction])

# Define trial and test functions
du = BlockTrialFunction(W)
ul = BlockFunction(W)
v = BlockTestFunction(W)
gap = Function(V2, name="Gap")

# Split trial and test functions into displacement and lagrange multiplier
(u, l) = block_split(ul)
(v_u, v_l) = block_split(v)

# Measure
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Boundary Conditions
# -----------------------------------------------------------------------------
u_bot = Expression(("0.0*x[0]", "0.0*x[1]", "0.0*x[2]"), degree=1)
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

# Fix bottom surface
bc_bot = DirichletBC(W.sub(0), u_bot, bot)
# Sub displacement subspace to obtain direction of interest
bc_r_l = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_r = DirichletBC(W.sub(0).sub(0), u_lr, sym_x)
bc_r_b = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)
bc_r_f = DirichletBC(W.sub(0).sub(2), u_bf, sym_z)

bc = BlockDirichletBC([[bc_bot, bc_r_l, bc_r_r, bc_r_b, bc_r_f], []])

# Constitutive relation
# -----------------------------------------------------------------------------
# Strain function
def eps(u):
    return sym(grad(u))
# Stress function defined using the Neo-Hookean Model
def sigma(u):
    return lmbda*tr(eps(u))*Identity(3) + 2.0*mu*eps(u)
# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.
# Define the augmented lagrangian
def aug_l(x):
    return l + k_pen*(u[1] - indent)

# Expression for indenter
indent = Expression(("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)"), \
                    depth=depth, R=R, degree=0)

# Weak form
F = [inner(sigma(u), eps(v_u))*dx + v_u[1]*ppos(aug_l(l))*ds,
    v_l*(u[1]-indent)*ds - (1/k_pen)*l*v_l*ds]
# Jacobian
J = block_derivative(F, ul, du)

# Solve with PETsc SNES solver
problem = BlockNonlinearProblem(F, ul, bc, J)
solver = BlockPETScSNESSolver(problem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

results = XDMFFile(name)
while steps <= tot_steps:
    print(steps)
    solver.solve()
    # Update
    steps += 1
    depth += 0.01
    indent.depth = depth
    # Split results
    (u, l) = ul.block_split()
    # Save results
    u.rename("Displacement", "u")
    l.rename("Lagrange Multiplier", "l")
    gap.assign(project(indent-u[1], V2))
    # Parameters will share the same mesh
    results.parameters["flush_output"] = True
    results.parameters["functions_share_mesh"] = True
    # Write
    results.write(u, steps)
    results.write(l, steps)
    results.write(gap, steps)