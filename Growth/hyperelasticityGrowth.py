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
userpar.add("de", 1)
userpar.add("nu", 0.50)
userpar.add("tot_steps", 500)
userpar.add("gamma", 1.1)
userpar.add("E_1", 50)
userpar.add("thick", 1)
userpar.add("w", 4/0.33*2*pi)
userpar.parse()

# Global stepping and displacement stepping parameters
steps = 1                           # Displacement step counter (updated within loop)
tot_steps = userpar["tot_steps"]    # Total displacement steps

# Initialize tolerance
tol = 1E-14
# Thickness and width
thick = userpar["thick"]
w = userpar["w"]

# Elasticity Parameters
E_1 = userpar["E_1"]
E_2 = Constant(1.0)

# Assign Young's Modulus to each domain
E = Expression('x[1] >= 100 - thick + tol ? E_1 : E_2', degree=0,
               tol=tol, thick=thick, E_1=E_1, E_2=E_2)

# Assign theta to each domain
# Expression for time step for updating in loop
# theta = Expression("growth", growth=growth, degree=0)
# growth factor
gamma = userpar["gamma"]        # Total growth factor
inc = (gamma-1)/tot_steps       # Incremental increase for growth per step

growth = 1
no_growth = 1
theta = Expression('x[1] >= 100 - thick + tol ? growth : no_growth',degree=0,
                    tol=tol, thick=thick, growth=growth, no_growth=no_growth)

# Solve for mu and lambda
nu = userpar["nu"]                    # Poisson Ratio
# lmbda = E*nu/((1 + nu)*(1 - 2*nu))  # Bulk Modulus
mu = E/(2*(1 + nu))                   # Shear Modulus
# alpha = nu/(1-2*nu)                 # Ratio of lmbda/2*mu

# Define Dirichlet boundary expressions
c = Constant((0.0))

# Define body force (B) and traction (T) terms
B  = Constant((0.0, 0.0))       # Body force per unit volume
T  = Constant((0.0, 0.0))       # Traction force on the boundary

# Create mesh and define function space
mesh = Mesh("hyperelasticity2D_gmsh_t_" + str(thick) + ".xml")

# Taylor-Hood space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Hydrostatic Pressure
TH = MixedElement([P2,P1])
V  = FunctionSpace(mesh, TH)

# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)

# The portions of the boundary on which Dirichlet boundary conditions will be
# applied are defined

# Mark boundary subdomains
# Subdomain 'left' corresponds to the part of the boundary on which 'x=0'
# left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
# # Subdomain 'right' corresponds to the part of the boundary on which 'x=1'
# right = CompiledSubDomain("near(x[0], side) && on_boundary", side = w)
# # Subdomain 'bottom' corresponds to the part of the boundary on which 'y=0'
# bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

class left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],0.0, 0.1)

class right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],w, 0.1)

class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1],0.0,0.1)

# Create class of all boundaries
class AllBoundaries(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Create boundaries
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Mark top boundary to restrict the lagrange multiplier
left = left()
right = right()
bottom = bottom()
allboundaries = AllBoundaries()

#
allboundaries.mark(sub_domains, 1)
left.mark(sub_domains, 0)
right.mark(sub_domains, 0)
bottom.mark(sub_domains, 0)

file_results = XDMFFile("SubDomains.xdmf")
file_results.write(sub_domains)

# Roller on left boundary
bcl = DirichletBC(V.sub(0).sub(0), c, left)
# Ramped displacement in x direction on right boundary
bcr = DirichletBC(V.sub(0).sub(0), c, right)
# Fixed at bottom boundary
bcb = DirichletBC(V.sub(0), Constant((0,0)), bottom)
# Combine boundary conditions
bcs = [bcb, bcl, bcr]

# Trial and test functions, and the most recent approximate displacement, u are
# defined on the finite element space V.

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
up  = Function(V)                 # Displacement

# NOTE that this is a shallow copy not a deep copy
(v_u, v_p) = split(v)            # Split the test functions
(u, p) = split(up)                # Split the solutions to (u,p)

# Define kinematic quantities involved in the model
# Kinematics
l = len(u)          # Length of displacement vector
I = Identity(l)     # Identity tensor

# Deformation Gradients
F = I + grad(u)
Fg = theta*I
Fe = F*inv(Fg)

# Cauchy-Green Tensor: C = F^T F
C = F.T*F           # Total
Ce = Fe.T*Fe        # Elastic

# Invariants of deformation tensors
Ic = tr(C) + 1      # Total
J  = det(F)         # Total
Ie = tr(Ce) + 1     # Elastic
Je = det(Fe)        # Elastic
Jg = det(Fg)        # Plastic - Growth

# Define terms for surface tension
ds = Measure("ds")(subdomain_data=sub_domains)   # Measure ds according to subdomains
N = FacetNormal(mesh)                            # Normal vector in the reference configuration
NansonOp = (cofac(F))                            # Element of area transformation operator
deformed_N = dot(NansonOp,N)                     # Surface element vector
Jsurf = sqrt(dot(deformed_N,deformed_N))         # Norm of the surface element vector

# First Piola-Kirchoff stress tensor
def P(u):
    return inv(Fg)*(J*mu*Fe-p*Je*inv(Fe.T))

# Total potential energy
F1 = (inner(P(u), grad(v_u)) + inner(Je-1., v_p))*dx(metadata={"quadrature_degree": 4}) \
     - dot(B, v_u)*dx - dot(T, v_u)*ds

# Variational problem where we have two equations for the weak form
de = userpar["de"]
# surface_energy_density = project(d, FunctionSpace(mesh,'DG',0))
# surface_energy = surface_energy_density*Jsurf*ds(1)
surface_energy_density = de*Jsurf
surface_energy = surface_energy_density*ds(1)
F2 = derivative(surface_energy, up, v)
WF = F1 + F2

# Compute Jacobian of F
Jacobian = derivative(WF, up, du)

problem = NonlinearVariationalProblem(WF, up, bcs, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file
name = "NeoHookean"
sim_param1 = "_w_%.2f" % (w)
sim_param2 = "_de_%.0f" % (de)
file_results = XDMFFile(name + sim_param1 + sim_param2 + ".xdmf")

# Solve for each value using the previous solution as a starting point
while (steps < tot_steps):
    # Print outs to track code progress
    print("Steps: " + str(steps))
    print("Time: " + str(growth))
    steps += 1                  # Update total steps
    # Update for Growth
    growth = growth+inc                 # Update growth with incremental growth value
    theta.growth = growth               # Update in expression class

    # Solve variational problem
    solver_problem.solve()

    # Note that this is now a deep copy not a shallow copy
    (u, p) = up.split()

    # Project nominal stress
    PTensor = project(P(u), TT)

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    PTensor.rename("Nominal Stress", "P")

    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True

    # Write to .xdmf results file
    file_results.write(u, steps)
    file_results.write(p, steps)
    file_results.write(PTensor,steps)
