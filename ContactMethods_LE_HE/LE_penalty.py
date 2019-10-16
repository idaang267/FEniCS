# 3D Linear elastic problem with penalty approach to contact

# Created to troubleshoot contact methods

# Import all modules
from dolfin import *
from mshr import *

# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": False}}

# Define classes for surfaces
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

# Convert classes
top = top()
bot = bot()
sym_x = symmetry_x()
sym_z = symmetry_z()

# Parameters
# -----------------------------------------------------------------------------
name = "LE_penalty.xdmf"
k_pen = 1e4                         # Penalty parameter
mesh_res = 20                       # Rough resolution of mesh
# Steps to control amount of indentation
steps = 0                           # Updated in loop
tot_steps = 10                      # Total amount of steps
# Define the indenter (obstacle) with a parabola
R = 0.25                            # Radius
depth = 0.00                        # Depth
# Material Parameters
E = Constant(10.0)                  # Young's Modulus
nu = Constant(0.3)                  # Poisson's Ratio
mu = E/(2*(1+nu))                   # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))      # Lamé Coefficient

# Create Mesh and define subdomains and boundaries of mesh
# -----------------------------------------------------------------------------
domain = Box(Point(0., 0., 0.), Point(1.0, 1.0, 1.0))
mesh = generate_mesh(domain, mesh_res)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
# Create boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Mark boundary to restrict ds, the surface, to the top surface
top.mark(boundaries, 1)

# Measure redefines the subdomains and boundaries
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Define all function spaces
# -----------------------------------------------------------------------------
V = VectorFunctionSpace(mesh, "Lagrange", 2)
# Function spaces for projection
V2 = FunctionSpace(mesh, "CG", 1)       # Gap
V0 = FunctionSpace(mesh, "DG", 0)       # Pressure

# Define trial and test functions and unknowns (displacement)
du = TrialFunction(V)
v_u = TestFunction(V)
u = Function(V)
# Define gap and pressure (name for paraview)
gap = Function(V2, name="Gap")
p = Function(V0, name="Contact pressure")

# Boundary Conditions
# -----------------------------------------------------------------------------
# Fixed surface
u_bot = Expression(("(0.0)*x[0]", "0.0*x[1]", "(0.0)*x[2]"), degree=1)
# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(0.0)*x[0]"), degree=1)
u_bf = Expression(("(0.0)*x[2]"), degree=1)

# Specify fixed bottom surface
bc_bot = DirichletBC(V, u_bot, bot)
# Sub displacement subspace to obtain direction of interest
bc_r_l = DirichletBC(V.sub(0), u_lr, sym_x)
bc_r_r = DirichletBC(V.sub(0), u_lr, sym_x)
bc_r_b = DirichletBC(V.sub(2), u_bf, sym_z)
bc_r_f = DirichletBC(V.sub(2), u_bf, sym_z)

# Combine all boundary conditions
bc = [bc_bot, bc_r_l, bc_r_r, bc_r_b, bc_r_f]

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

# Weak form: one equation
F = inner(sigma(u), eps(v_u))*dx + k_pen*v_u[1]*ppos(u[1]-indent)*ds(1)

# Jacobian - Directional derivative
J = derivative(F, u, du)

#  Setup Non-linear variational problem using the PETSc SNES solver
problem = NonlinearVariationalProblem(F, u, bc, J=J)
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file since we have multiple fields
results = XDMFFile(name)

# Solve with Newton solver using the previous solution as a starting point
while steps <= tot_steps:
    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Depth: " + str(depth))

    # Solve the nonlinear problem
    solver.solve()

    # Update all parameters
    steps += 1                  # Update steps
    depth += 0.01               # Update the depth of the indenter
    indent.depth = depth        # Update depth in expression class

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    gap.assign(project(indent-u[1], V2))
    p.assign(-project(sigma(u)[1, 1], V0))
    # Parameters will share the same mesh
    results.parameters["flush_output"] = True
    results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    results.write(u, steps)
    results.write(gap, steps)
    results.write(p, steps)
