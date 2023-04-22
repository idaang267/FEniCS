# Hyperelasticity
#==============================================================================
# This example demonstrates the solution of a three-dimensional elasticity
# problem.
# First, the required modules are imported
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
# UFLACS is a collection of algorithms for processing symbolic unified form
# languages forms and expressions
parameters["form_compiler"]["representation"] = "uflacs"

# All parameters
# Control the number of loop iterations through user parameters
# Parse parameters from command line
parameters.parse()
UserPar = Parameters("user")
UserPar.add("Model", "Incomp")
UserPar.add("TotSteps", 20)     # Total length of the simulation
UserPar.add("OutPlaneN", 15)
UserPar.add("InPlaneN", 15)
UserPar.add("nu", 0.49)
# Add user parameters in the global parameter set
UserPar.parse()

#------------------------------------------------------------------------------
Model = UserPar["Model"]                   # Keywords: Comp or Incomp (Penalty)
TotSteps = UserPar["TotSteps"]              # Total simulation steps
RampArray = np.linspace(0, 0.5, TotSteps)   # For expressions
ThetaArray = np.linspace(0, pi/3, TotSteps) # For expressions
# Elasticity parameters
E  = 10.0                                   # Young's Modulus
nu = UserPar["nu"]                          # Poisson's ratio
mu = Constant(E/(2*(1 + nu)))               # Lame parameter
# Note that lambda is a reserved keyword in Python, hence the misspelling 'lmbda'
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu))) # Lame parameter
kappa = Constant(E/(3*(1 - 2*nu)))           # Bulk Modulus
# Two constants are declared for the body force (B) and traction (T) terms
B  = Constant((0.0, 0.0, 0.0))         # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))         # Traction force on the boundary Gamma_N
# Define Dirichlet boundary (x = 0 or x = 1) expressions
c = Constant((0.0, 0.0, 0.0))           # More efficent to use Constant not Expression
r = Expression(("scale*0.0",
                "scale*(0.5 + (x[1] - 0.5)*cos(theta) - (x[2] - 0.5)*sin(theta) - x[1])",
                "scale*(0.5 + (x[1] - 0.5)*sin(theta) + (x[2] - 0.5)*cos(theta) - x[2])"),
                scale = 0, theta = 0, degree=2)
# Number of elements in-plane and out-plane
InPlaneN, OutPlaneN = UserPar["InPlaneN"], UserPar["OutPlaneN"]
SimPar = "S_%.0d_Nx_%.0d_N_%.0d_nu_%.1e" % (TotSteps, OutPlaneN, InPlaneN, nu)
print(Model)
if Model == "Comp":
    SaveDir = '../Results/EnergyComp/'
else:
    SaveDir = '../Results/EnergyPenalty/'
File(SaveDir+ SimPar +"/parameters.xml") << UserPar

# Create tetrahedral mesh of the domain and a function space on this mesh
mesh = UnitCubeMesh(OutPlaneN, InPlaneN, InPlaneN) # Domain - Omega

# On this mesh, we define a function space of continuous piecewise linear
# vector polynomials
V = VectorFunctionSpace(mesh, "Lagrange", 1)
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)
DG0 = FunctionSpace(mesh, "DG", 0)

# Trial and test functions, and the most recent approximate displacement, u are
# defined on the finite element space V.
# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V, name='Displacement') # Displacement from previous iteration

# Kinematics
#------------------------------------------------------------------------------
d = len(u)          # Length of displacement vector
I = Identity(d)     # Identity tensor
F = I + grad(u)     # Deformation gradient
C = F.T*F           # Right Cauchy-Green tensor C = F^T F

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# The portions of the boundary on which Dirichlet boundary conditions will be
# applied are defined

# Gamma_{D_0}: The boundary subdomain 'left' corresponds to the part of the
# boundary on which 'x=0'
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
# Gamma_{D_1}: The boundary subdomain 'right' corresponds to the part of the
# boundary on which 'x=1'
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# The boundary subdomains and the boundary condition expressions are collected
# together in two Dirichlet BC objects, one for each part of the boundary
bcl = DirichletBC(V, c, left)       # Gamma_{D_0}
bcr = DirichletBC(V, r, right)      # Gamma_{D_1}
bcs = [bcl, bcr]                    # Combine all boundary conditions

# The strain energy density and the total potential energy are defined (using UFL syntax)
#------------------------------------------------------------------------------
# Stored strain energy density (compressible neo-Hookean model)
def psi(Ic, J):
    if Model == "Comp":
        return (mu/2)*(Ic - 3-  mu*ln(J)) + (lmbda/2)*(ln(J))**2
    else:
        return (mu/2)*(Ic - 3 - 2*ln(J)) + (kappa/2)*(J-1)**2

# 1st Piola-Kirchoff Stress
def P(F):
    if Model == "Comp":
        return mu*(F-inv(F.T)) + lmbda*ln(J)*inv(F.T)
    else:
        return mu*(F-inv(F.T)) + kappa*(J-1)*J*inv(F.T)

# Total potential energy
Pi = psi(Ic,J)*dx - dot(B, u)*dx - dot(T, u)*ds

# Directional derivatives are now computed of Pi and L
# Compute first variation of Pi, directional derivative about u in the direction of v
WF = derivative(Pi, u, v)
Jc = derivative(WF, u, du)        # Compute Jacobian of F

# Finally, the solution 'u' is saved to a file name
file = XDMFFile(SaveDir + SimPar + "/Results.xdmf"); # Save solution in VTK format
file.parameters["flush_output"] = True
file.parameters["functions_share_mesh"] = True

# Post processing
PostProc = np.zeros((TotSteps, 5))

for (Step, Scale) in enumerate(RampArray):
    # Print outs to track code progress
    print("\033[1;32m--- Step {0:2d}: scale = {1:2f} ---\033[1;m".format(Step, Scale))

    # The complete variational problem can now be solved by a single call to solve
    solve(WF == 0, u, bcs, J=Jc)
    # Project nominal stress to Discontinuous Galerkin space
    PK = project(P(F), T_DG0)
    DetF = project(J, DG0)
    # Rename files for visualization
    PK.rename("Piola Kirchoff Stress", "PK")
    DetF.rename("J", "DetF")
    # Save all parameters of interest
    file.write(u, Step)
    file.write(PK, Step)
    file.write(DetF, Step)
    # Postprocessing
    PostProc[Step] = np.array([Step, Scale, u(1,1,1)[1], u(1,1,1)[2], DetF(1,1,1)])
    np.savetxt(SaveDir + SimPar + '/PostProc.txt', PostProc)

    # Update expressions
    r.scale = Scale
    r.theta = ThetaArray[Step]
