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
name = "LE_penalty.xdmf"
# Steps
steps = 0
tot_steps = 10
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

# Create Mesh
domain = Box(Point(0., 0., 0.), Point(1.0, 1.0, 1.0))
mesh = generate_mesh(domain, 20)

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

# Function space
V = VectorFunctionSpace(mesh, "Lagrange", 2)

# Define trial and test functions
du = TrialFunction(V)
u = Function(V)
v_u = TestFunction(V)

# Measure
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Boundary Conditions
# Lateral surfaces
# left and right surfaces
class symmetry_x(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) or near(x[0], 1.0) and on_boundary
# front and back surfaces
class symmetry_z(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0.0) or near(x[2], 1.0) and on_boundary

sym_x = symmetry_x()
sym_z = symmetry_z()

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(0.0)*x[0]"), l0=l0, degree=1)
u_bf = Expression(("(0.0)*x[2]"), l0=l0, degree=1)

u_bot = Expression(("(0.0)*x[0]", "0.0*x[1]", "(0.0)*x[2]"), l0=l0, degree=1)
bc = DirichletBC(V, u_bot, bot)

# # Initial Conditions (IC)
# #------------------------------------------------------------------------------
# # Initial conditions are created by using the class defined and then
# # interpolating into a finite element space
# block_u = ul.sub(0)
# block_l = ul.sub(1)
# # Interpolate on subfunction
# u_ini = Expression(("(l0-1)*x[0]", "(l0-1)*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
# block_u.interpolate(u_ini)
# block_l.interpolate(Constant(0.0))
# # Update the original block functions by assign(receiving_funcs, assigning_func)
# # Current solution
# assign(ul.sub(0), block_u)
# assign(ul.sub(1), block_l)

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

# Indenter
indent = Expression(("-depth+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*R)"), \
                    depth=depth, R=R, degree=0)
# # Indenter
# indent = Expression(("-depth"), depth=depth, degree=0)
# Weak form
F = inner(sigma(u), eps(v_u))*dx + 1E4*v_u[1]*ppos(u[1]-indent)*ds(1)

# Jacobian
J = derivative(F, u, du)

# SNES solver > Setup Non-linear variational problem
problem = NonlinearVariationalProblem(F, u, bc, J=J)
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)

results = XDMFFile(name)
while steps <= tot_steps:
    print(steps)
    solver.solve()
    # Update
    steps += 1
    depth += 0.01
    indent.depth = depth
    # Save results
    u.rename("Displacement", "u")
    # Parameters will share the same mesh
    results.parameters["flush_output"] = True
    results.parameters["functions_share_mesh"] = True
    # Write
    results.write(u, steps)
