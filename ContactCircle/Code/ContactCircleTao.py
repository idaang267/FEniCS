# Contact Mechanics Problem using PETSc's TAO solver for nonlinear
# (bound-constrained) optimisation problems

# This example considers a heavy hyperelastic circle under body force in a box
# of the same size

# Copyright (C) 2012 Corrado Maurini
# Modified by Tianyi Li 2014
# Modified by Ida Ang 2023

# Import Modules
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

# Check if PETSc is installed otherwise exit
if not has_petsc():
    print("DOLFIN must be compiled with PETSc to run this demo.")
    exit(0)

# Define the minimization problem by using the OptimisationProblem class
class ContactProblem(OptimisationProblem):
    def __init__(self):
        OptimisationProblem.__init__(self)
    # Objective function
    def f(self, x):
        u.vector()[:] = x
        return assemble(Pi)
    # Gradient of the objective function
    def F(self, b, x):
        u.vector()[:] = x
        assemble(F, tensor=b)
    # Hessian of the objective function
    def J(self, A, x):
        u.vector()[:] = x
        assemble(J, tensor=A)

# Create the PETScTAOSolver
solver = PETScTAOSolver()

# Set parameters for the solver
solver.parameters["method"] = "tron"
# solver.parameters["linear_solver"] = "nash"
solver.parameters["line_search"] = "gpcg"
    # when using gpcg make sure that you have a constant Hessian
# solver.parameters["preconditioner"] = "ml_amg"
solver.parameters["monitor_convergence"] = True
solver.parameters["report"] = False
# Uncomment this line to see the available parameters
#info(parameters, True)

# Symmetry condition (to block rigid body rotations)
def symmetry_line(x, on_boundary):
    return near(x[0], 0)

# Parse (PETSc) parameters
parameters.parse()

# Elasticity Parameters
E = 10.0                                        # Young's Modulus
nu = 0.3                                        # Poisson Ratio
mu = Constant(E/(2*(1 + nu)))                   # Lamé Parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))    # Lamé Parameter

# Body force per unit volume, updated as we go
B = Expression(("0.0", "BodyInc"), BodyInc = 0.0, degree = 0 )
Body = np.linspace(0, -1, 100)

# Read pre-made mesh
Dir = "../Results/TAO/"
mesh = Mesh("../Mesh/CircleXY.xml")

# Create function space
V = VectorFunctionSpace(mesh, "Lagrange", 1)
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)      # Tensor function space for projection

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

# Kinematics
d = len(u)                       # Spatial dimension
I = Identity(d)                  # Identity tensor
Fdef = I + grad(u)                  # Deformation gradient
C = Fdef.T*Fdef                        # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
Jdef  = det(Fdef)

# The displacement u must be such that the current configuration
# doesn't escape the constraint of a box
#-------------------------------------------------------------------------------
# Subclass funciton space to constrain x (0) direction
bc = DirichletBC(V.sub(0), Constant(0.0), symmetry_line)
# Identifies the right side and top sides of the square
constraint_l = Expression(("-1.-x[0]", "x[0] < -0.2 ? -1. -x[1] : x[0] > 0.2 ? -1-x[1] : -1-x[1]-inc"), \
                           inc=5.0,  degree=1)
# Identifies the left side and bottom side of the square
constraint_u = Expression(("1.-x[0]", "1.-x[1]"), degree=1)

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic-2 - 2*ln(Jdef)) + (lmbda/2)*(ln(Jdef))**2

def P(u):
    return mu*(Fdef - inv(Fdef.T)) + lmbda*(ln(Jdef))*Jdef*inv(Fdef.T)

# Total potential energy
Pi = psi*dx - dot(B,u)*dx

# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)
# Compute Jacobian of F
J = derivative(F, u, du)

# Tao solver minimum and maximum
u_min = interpolate(constraint_l, V)
u_max = interpolate(constraint_u, V)

# BC will be incorporated into the lower (min) and upper (max) bounds
bc.apply(u_min.vector())
bc.apply(u_max.vector())

# Save solution in XDMF format if available
file = XDMFFile(mesh.mpi_comm(), Dir + "Results.xdmf")
file.parameters["rewrite_function_mesh"] = False
file.parameters["functions_share_mesh"]  = True     # Allows saving of multiple functions
file.parameters["flush_output"]          = True     # Saves the file in case of interruption

# Looping with the update of body force
for (i, inc) in enumerate(Body):
    print("Increment of Body Force: ", Body[i])

    # Solve the problem for the specific step
    solver.solve(ContactProblem(), u.vector(), u_min.vector(), u_max.vector())

    # Project nominal stress function to tensor space
    PTensor = project(P(u), T_DG0)

    # Rename for visualization in paraview
    u.rename("Disp", "u")
    PTensor.rename("Nominal Stress", "P")
    # Write to xdmf file
    file.write(u, i)
    file.write(PTensor, i)

    # Update body force before next iteration
    B.BodyInc = Body[i]
