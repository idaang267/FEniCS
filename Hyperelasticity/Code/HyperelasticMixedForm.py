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
# parameters["form_compiler"]["log_level"] = INFO
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

# Parse parameters from command line
parameters.parse()
# Control the number of loop iterations through user parameters
UserPar = Parameters("user")
UserPar.add("TotSteps", 20)     # Displacement step number
UserPar.add("OutPlaneN", 10)
UserPar.add("InPlaneN", 10)
UserPar.add("nu", 0.4)
UserPar.parse()

# All parameters
#------------------------------------------------------------------------------
TotSteps = UserPar["TotSteps"]
RampArray = np.linspace(0, 0.5, TotSteps)
ThetaArray = np.linspace(0, pi/3, TotSteps)
# Elasticity parameters
E = 10.0                          # Young's Modulus, Poisson's Ratio
nu = UserPar["nu"]                          # Poisson's ratio
mu = Constant(E/(2*(1 + nu)))                # Lame Parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu))) # Lame parameter
# Definitions of constants from hyperelasticity demo
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
# Define Dirichlet boundary (x = 0 or x = 1) expressions
c = Constant((0.0, 0.0, 0.0))
r = Expression(("scale*0.0",
                "scale*(.5 + (x[1] - .5)*cos(theta) - (x[2] - .5)*sin(theta) - x[1])",
                "scale*(.5 + (x[1] - .5)*sin(theta) + (x[2] - .5)*cos(theta) - x[2])"),
                scale=0, theta=0, degree=2)
# Number of elements in-plane and out-plane
InPlaneN, OutPlaneN = UserPar["InPlaneN"], UserPar["OutPlaneN"]
SimPar = "S_%.0d_Nx_%.0d_N_%.0d_nu_%.1e" % (TotSteps, OutPlaneN, InPlaneN, nu)
SaveDir = '../Results/WeakFormUp/'

# Define mesh and mixed function space
mesh = UnitCubeMesh(OutPlaneN, InPlaneN, InPlaneN)
# Taylor-Hood space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Hydrostatic Pressure
TH = MixedElement([P2,P1])
V  = FunctionSpace(mesh, TH)
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)
DG0 = FunctionSpace(mesh, "DG", 0)

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, hydrostatic pressure)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
w  = Function(V)                 # Solutions (u,p) from previous iteration

# NOTE that this is a shallow copy not a deep copy
(v_u, v_p) = split(v)            # Split the test functions
(u, p) = split(w)                # Split the solutions to (u,p)

# Kinematics
#------------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark boundary subdomains
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define Dirichlet boundary (x = 0 or x = 1)
# The Dirichlet BCs are specified in a subspace
bcl = DirichletBC(V.sub(0), c, left)
bcr = DirichletBC(V.sub(0), r, right)
bcs = [bcl, bcr]

# Where P = dW/dF:
def P(F):
    return mu*F + p*J*inv(F.T)

# Specify the quadrature degree for efficiency
WF = (inner(P(F), grad(v_u)) + inner(J-1., v_p))*dx(metadata={"quadrature_degree": 4}) \
    - dot(B, v_u)*dx - dot(T, v_u)*ds

# Directional derivative
J_o = derivative(WF, w, du)

# Solve variational problem
varproblem = NonlinearVariationalProblem(WF, w, bcs, J=J_o)
solver = NonlinearVariationalSolver(varproblem)
solver.parameters.update(snes_solver_parameters)
info(solver.parameters, False)

# Save results to an .xdmf file since we have multiple fields
file = XDMFFile(SaveDir + SimPar + "/Results.xdmf")
file.parameters["flush_output"] = True
file.parameters["functions_share_mesh"] = True

# Post processing
PostProc = np.zeros((TotSteps, 5))

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
for (Step, Scale) in enumerate(RampArray):
    # Print outs to track code progress
    print("\033[1;32m--- Step {0:2d}: scale = {1:2f} ---\033[1;m".format(Step, Scale))

    # Solve the nonlinear problem (using Newton solver)
    solver.solve()

    # Note that this is now a deep copy not a shallow copy
    (u, p) = w.split()

    # Project nominal stress to Discontinuous Galerkin space
    PK = project(P(F), T_DG0)
    DetF = project(J, DG0)

    PK.rename("Piola Kirchoff Stress", "PK")
    u.rename("Displacement", "u")
    p.rename("Pressure","p")
    DetF.rename("J", "DetF")

    file.write(u, Step)
    file.write(p, Step)
    file.write(PK, Step)
    file.write(DetF, Step)

    PostProc[Step] = np.array([Step, Scale, u(1,1,1)[1], u(1,1,1)[2], DetF(1,1,1)])
    np.savetxt(SaveDir + SimPar + '/PostProc.txt', PostProc)

    # Update the displacement value
    r.scale = Scale
    r.theta = ThetaArray[Step]
