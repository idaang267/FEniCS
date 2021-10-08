# Hyperelasticity
# Adapted from Hyperelasticity Demo:

# Using the Incompressible Hyperelasticity Strain Energy Formulation
#   The stored strain energy density:
#   Compressible neo-Hookean model:
#       W = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2
#   Incompressible neo-Hookean model (J = det(F) = 1):
#       W = (mu/2)*(Ic-3)
#   Add a lagrange multiplier term to enforce incompressibility:
#       W = (mu/2)*(Ic - 3) + p*(J-1)

# Edited to use SNES solver

# Import modules
from dolfin import *
import matplotlib.pyplot as plt     # For visualization
import numpy as np

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"
# The following two fields were in the original hyperelasticity demo
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
    # For processing UFLACS - unified form language expressions
parameters["form_compiler"]["log_level"] = INFO
    # Show progress of compiler
parameters["allow_extrapolation"] = True
# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          'absolute_tolerance':1e-5,
                                          'relative_tolerance':1e-5,
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

solver_par = NonlinearVariationalSolver.default_parameters()
solver_par.rename("solver")

# Control the number of loop iterations through user parameters
user_par = Parameters("user")
user_par.add("u_min",0.)        # Displacement minimum
user_par.add("u_max",0.5)       # Displacement maximum
user_par.add("u_nsteps",5)     # Displacement step number

# Add user parameters in the global parameter set
parameters.add(user_par)

# Parse parameters from command line
parameters.parse()
info(parameters,True)
user_par = parameters.user

# Define mesh and mixed function space
mesh = UnitCubeMesh(20, 10, 10)
# Taylor-Hood space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Hydrostatic Pressure
TH = MixedElement([P2,P1])
V  = FunctionSpace(mesh, TH)

# Mark boundary subdomains
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define Dirichlet boundary (x = 0 or x = 1)
c = Constant((0.0, 0.0, 0.0))
r = Expression(("scale*0.0",
                "scale*(scale + (x[1] - scale)*cos(theta) - (x[2] - scale)*sin(theta) - x[1])",
                "scale*(scale + (x[1] - scale)*sin(theta) + (x[2] - scale)*cos(theta) - x[2])"),
                scale = 0, theta = pi/15, degree=2)
# The Dirichlet BCs are specified in a subspace
bcl = DirichletBC(V.sub(0), c, left)
bcr = DirichletBC(V.sub(0), r, right)
bcs = [bcl, bcr]

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, hydrostatic pressure)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
w  = Function(V)                 # Solutions (u,p) from previous iteration

# NOTE that this is a shallow copy not a deep copy
(v_u, v_p) = split(v)            # Split the test functions
(u, p) = split(w)                # Split the solutions to (u,p)

# Definitions of constants from hyperelasticity demo
B  = Constant((0.0, -0.5, 0.0))  # Body force per unit volume
T  = Constant((0.1,  0.0, 0.0))  # Traction force on the boundary

# Kinematics
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
E, nu = 10.0, 0.3
mu = Constant(E/(2*(1 + nu)))

# Where P = dW/dF:
def P(u):
    return mu*F + p*J*inv(F.T)

# Specify the quadrature degree for efficiency
F = (inner(P(u), grad(v_u)) + inner(J-1., v_p))*dx(metadata={"quadrature_degree": 4}) \
    - dot(B, v_u)*dx - dot(T, v_u)*ds

# Directional derivative
J_o = derivative(F, w, du)

# Solve variational problem
varproblem = NonlinearVariationalProblem(F, w, bcs, J=J_o)
solver = NonlinearVariationalSolver(varproblem)
solver.parameters.update(snes_solver_parameters)
info(solver.parameters, True)
(iter, converged) = solver.solve()

# Loading parameter (list of values of surface tension for the simulations)
u_list = np.linspace(user_par.u_min,user_par.u_max,user_par.u_nsteps)

# Save results to an .xdmf file since we have multiple fields
file_results = XDMFFile("results.xdmf")

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
for scale in u_list:
    # Update the displacement value
    r.scale = scale
    # Solve the nonlinear problem (using Newton solver)
    solver.solve()

    # Note that this is now a deep copy not a shallow copy
    (u, p) = w.split()
    u.rename("Displacement", "u")
    p.rename("Pressure","p")
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    file_results.write(u,scale)
    file_results.write(p,scale)
