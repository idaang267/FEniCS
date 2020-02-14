# Hyperelasticity
# ===============
# Background
# ----------
# This example demonstrates the solution of a two-dimensional elasticity
# problem.
#
# Equation and problem definition
# -------------------------------

# First, the required modules are imported
from dolfin import *
import matplotlib.pyplot as plt
from ufl import cofac, rank

# Optimization options for the form compiler
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 50,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "error_on_nonconvergence": False}}

# Set the user parameters, can be parsed from command line
parameters.parse()
userpar = Parameters("user")
userpar.add("gamma", 0.05)
userpar.add("nu", 0.3)
userpar.add("tot_steps", 19)
userpar.parse()

# Global stepping and displacement stepping parameters
steps = 1                           # Displacement step counter (updated within loop)
tot_steps = userpar["tot_steps"]    # Total displacement steps

# Time parameters
dt = 1                       # Starting time step
# Expression for time step for updating in loop
theta = Expression("dt", dt=dt, degree=0)
t = 0.0                         # Initial time (updated in loop)
c_exp = 1.05                     # Control the time step increase
# Initialize tolerance
tol = 1E-14
# Surface Energy: Gamma term
gamma = userpar["gamma"]
Gamma = Expression("gamma", gamma=gamma, degree=0)

# Elasticity Parameters
E_1 = Constant(100.0)
E_2 = Constant(1.0)
# Assign Young's Modulus to each domain
E = Expression('x[1] >= 0.09 + tol ? E_1 : E_2', degree=0,
               tol=tol, E_1=E_1, E_2=E_2)

# Solve for mu and lambda
nu = userpar["nu"]                  # Poisson Ratio
lmbda = E*nu/((1 + nu)*(1 - 2*nu))  # Bulk Modulus
mu = E/(2*(1 + nu))                 # Shear Modulus
alpha = nu/(1-2*nu)                 # Ratio of lmbda/2*mu

# Define Dirichlet boundary expressions
c = Constant((0.0))

# Define body force (B) and traction (T) terms
B  = Constant((0.0, 0.0))       # Body force per unit volume
T  = Constant((0.0, 0.0))       # Traction force on the boundary

name = "BlatzKo"
sim_param1 = "_g_%.2f" % (gamma)
sim_param2 = "_nu_%.2f" % (nu)

# Create mesh and define function space
mesh = Mesh("hyperelasticity2D_gmsh.xml")

# On this mesh, we define a function space of continuous piecewise linear
# vector polynomials (a Lagrange vector element space)
V = VectorFunctionSpace(mesh, "Lagrange",1)
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)

# The portions of the boundary on which Dirichlet boundary conditions will be
# applied are defined

# Mark boundary subdomains
# Subdomain 'left' corresponds to the part of the boundary on which 'x=0'
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
# Subdomain 'right' corresponds to the part of the boundary on which 'x=1'
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
# Subdomain 'bottom' corresponds to the part of the boundary on which 'y=0'
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

# Create class of all boundaries
class AllBoundaries(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Create boundaries
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Mark top boundary to restrict the lagrange multiplier
allboundaries = AllBoundaries()
allboundaries.mark(sub_domains, 1)
left.mark(sub_domains, 0)
right.mark(sub_domains, 0)
bottom.mark(sub_domains, 0)

# Roller on left boundary
bcl = DirichletBC(V.sub(0), c, left)
# Ramped displacement in x direction on right boundary
bcr = DirichletBC(V.sub(0), c, right)
# Roller on bottom boundary
bcb = DirichletBC(V.sub(1), c, bottom)
# Combine boundary conditions
bcs = [bcl, bcr, bcb]

# Trial and test functions, and the most recent approximate displacement, u are
# defined on the finite element space V.

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement

# Define kinematic quantities involved in the model
# Kinematics
d = len(u)          # Length of displacement vector
I = Identity(d)     # Identity tensor

# Deformation Gradients
F = I + grad(u)
Fg = theta*I
Fe = F*inv(Fg)

# Cauchy-Green Tensor: C = F^T F
C = F.T*F           # Total
Ce = Fe.T*Fe        # Elastic

# Invariants of deformation tensors
Ic = tr(C) + 1     # Total
J  = det(F)        # Total
Ie = tr(Ce) + 1    # Elastic
Je = det(Fe)       # Elastic
Jg = det(Fg)       # Plastic - Growth

# Define terms for surface tension
ds = Measure("ds")(subdomain_data=sub_domains)   # Measure ds according to subdomains
N = FacetNormal(mesh)                            # Normal vector in the reference configuration
NansonOp = (cofac(F))                            # Element of area transformation operator
deformed_N = dot(NansonOp,N)                     # Surface element vector
Jsurf = sqrt(dot(deformed_N,deformed_N))         # Norm of the surface element vector

# First Piola-Kirchoff stress tensor
def P(u):
    return mu*Fe - mu/2*Je**(-alpha)*inv(Fe.T)
#    return mu*Jg*(Fe - 0.5*Je**(-alpha)*inv(Fe.T))*inv(Fg.T)
#    return mu*Fe + (lmbda*ln(Je) - mu)*inv(Fe.T)

# Total potential energy
F1 = inner(P(u), grad(v))*dx - dot(B, v)*dx - dot(T, v)*ds

# Variational problem where we have two equations for the weak form
surface_energy_density = Gamma*Jsurf
surface_energy = surface_energy_density*ds(1)
F2 = derivative(surface_energy, u, v)
WF = F1 + F2

# Compute Jacobian of F
Jacobian = derivative(WF, u, du)

problem = NonlinearVariationalProblem(WF, u, bcs, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file
file_results = XDMFFile(name + sim_param1 + sim_param2 + ".xdmf")

# Solve for each value using the previous solution as a starting point
while (steps < tot_steps):
    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Time: " + str(dt))

    steps += 1                  # Update total steps
    dt = dt*c_exp               # Update time step with exponent value
    theta.dt = dt               # Update steps in expression class

    # Solve variational problem
    solver_problem.solve()

    # Project nominal stress
    PTensor = project(P(u), TT)

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    PTensor.rename("Nominal Stress", "P")

    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True

    # Write to .xdmf results file
    file_results.write(u, steps)
    file_results.write(PTensor,steps)
