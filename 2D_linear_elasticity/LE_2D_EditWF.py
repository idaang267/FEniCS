# Small strain solution for a 2D isotropic linear elastic medium in plane
# strain, in a traditional displacment-based finite element formulation

from dolfin import *
from ufl import rank

# Material Parameters
E = Constant(1e5)               # Young's Modulus
nu = Constant(0.3)              # Poisson's Ratio
mu = E/(2*(1+nu))               # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))  # Lamé Coefficient
# Body force
rho_g = 10
f = Constant((0,-rho_g))
# Traction
t = Constant((0.0,0.0))
# steps
steps = 0
tot_steps = 10
# Initial indenter depth of indentation
depth = 0.05
pen = Constant(1e4)

# Define mesh
# ---------------------
# Cantiliver beam modeled as a 2D medium of dimensions (L x H)
L = 1.0
H = 1.0
# Mesh density parameters
Nx = 10
Ny = 10
# Mesh density defined with a criss-crossed structured mesh
mesh = RectangleMesh(Point(0, 0), Point(L, H), Nx, Ny, "crossed")
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
top = CompiledSubDomain("near(x[1], 1.0) && on_boundary")
# Marking for visualization
boundaries.set_all(0)
top.mark(boundaries, 1)
# Measure redefines ds
ds = Measure('ds', subdomain_data=boundaries)

#XDMFFile("boundaries.xdmf").write(boundaries)

# Define the function space
V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
Vsig = TensorFunctionSpace(mesh, "DG", degree=0)

# Boundary Condition
# ---------------------
bot = CompiledSubDomain("near(x[1], 0.0) && on_boundary")
indent = Expression(("-depth"), depth=depth, degree=1)

bc_bot = DirichletBC(V, Constant((0.,0.)), bot)
bc_top = DirichletBC(V.sub(1), indent, top)
bc = [bc_bot, bc_top]

# Constitutive relation
# ---------------------
# Strain operator computes the 2x2 plane components of the symmetrized gradient
# tensor of any 2D vectorial field
def eps(v):
    return sym(grad(v))
# Stress operator computed
def sigma(v):
    return lmbda*Identity(2)*tr(eps(v)) + 2.0*mu*eps(v)
# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.

# Variational Problem
# ---------------------
du = TrialFunction(V)                   # Trial Function
v = TestFunction(V)                     # Test Function
u = Function(V, name="Displacement")    # Function named for paraview output
sig = Function(Vsig, name="Stress")

form = inner(sigma(u),eps(v))*dx - inner(f,v)*dx - inner(t,v)*ds \
        - pen*inner(v[1],u[1]-indent)*ds(1)

# Calculate Jacobian
J = derivative(form, u, du)

problem = NonlinearVariationalProblem(form, u, bc, J=J)
solver = NonlinearVariationalSolver(problem)

# Fields can be exported in a suitable format for visualization in Paraview
file_results = XDMFFile("elasticity_results.xdmf")

while (steps <= tot_steps):

    solver.solve()

    # Stress
    sig.assign(project(sigma(u), Vsig))
    # Control certain parameters
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Only one time step, 0
    file_results.write(u, steps)
    file_results.write(sig, steps)
    # Update total steps
    steps += 1
    depth += 0.05
    indent.depth = depth
