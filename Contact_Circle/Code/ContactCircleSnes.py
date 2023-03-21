# Contact mechanics problem using SNES solver for variational inequalities.

# This example considers a heavy hyperelastic circle under body force in a box
# of the same size

# Copyright (C) 2012 Corrado Maurini
# Modified by Corrado Maurini 2013

# Import Modules
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

# Check if PETSc is installed otherwise exit
if not has_petsc():
    print("DOLFIN must be compiled with PETSc to run this demo.")
    exit(0)

# Define the solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          "maximum_iterations": 50,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

# Subdomain for full boundary
class AllBound(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class MidPoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0, 0.01) and near(x[1], 0, 0.01)

# Symmetry condition (to block rigid body rotations)
class SymLine(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0]) < DOLFIN_EPS

AllBound = AllBound()
MidPoint = MidPoint()
SymLine = SymLine()

# Elasticity parameters
E = 10.0                                        # Young's Modulus
nu = 0.3                                        # Poisson Ratio
mu = Constant(E/(2*(1 + nu)))                   # Lamé Parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))    # Lamé Parameter

# Body force per unit volume, updated as we go
B = Expression(("0.0", "BodyInc"), BodyInc = 0.0, degree = 0 )
Body = np.linspace(0, -1, 50)

# Read pre-made mesh with format .xml
mesh = Mesh("Mesh/CircleXY.xml")

# Marking subdomains, lines or points
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
lines      = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
points     = MeshFunction("size_t", mesh, mesh.topology().dim()-2)

# Mark for visualization in Paraview
lines.set_all(0)
AllBound.mark(lines, 1)
file = XDMFFile("../Results/Lines.xdmf")
file.write(lines)

points.set_all(0)
MidPoint.mark(points, 1)
file = XDMFFile("../Results/Points.xdmf")
file.write(points)

# Create function spaces
V = VectorFunctionSpace(mesh, "Lagrange", 1)
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)      # Tensor function space for projection

# Define trial, test, and unknown functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement

# Kinematics
d = len(u)                      # Spatial Dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# The displacement u must be such that the current configuration (x+u)
# remains in the constraints
#-------------------------------------------------------------------------------s
# Subclass funciton space to constrain x(0) direction
bc = DirichletBC(V.sub(0), 0.0, SymLine, method="pointwise")

# Identifies the right side and top sides of the square
constraint_u = Expression(("1.-x[0]", "1.- x[1]"), tol=DOLFIN_EPS, degree=1)
# Identifies the left side and bottom side of the square
constraint_l = Expression(("-1.-x[0]", \
                           "x[0] < -0.2 ? -1. -x[1] : x[0] > 0.2 ? -1-x[1] : -1-x[1]-inc"), \
                           inc=5.0, degree=1)

# Identifies the left side and bottom side of the square
# constraint_l = Expression(("xmin - x[0]","ymin - x[1]"), \
#                             xmin=-1.0-DOLFIN_EPS, ymin=-1.0, degree=1)
# # Identifies the right side and top sides of the square

# Stored strain energy density (compressible neo-Hookean model)
def psi(u):
    return (mu/2)*(Ic - 2 - 2*ln(J)) + (lmbda/2)*(ln(J))**2

def P(u):
    return mu*(F - inv(F.T)) + lmbda*(ln(J))*inv(F.T)

# Total potential energy
Pi = psi(u)*dx - dot(B, u)*dx

# Compute first variation of Pi (directional derivative about u in the direction of v)
Fu = derivative(Pi, u, v)
# Compute Jacobian of F
Ju = derivative(Fu, u, du)

# SNES solver minimum and maximum
umin = interpolate(constraint_l, V)
umax = interpolate(constraint_u, V)

# Set up the non-linear problem
problem = NonlinearVariationalProblem(Fu, u, bc, J=Ju)
problem.set_bounds(umin, umax)

# Set up the non-linear solver
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)
info(solver.parameters, False)

# Save solution in XDMF format
file = XDMFFile("../Results/Results.xdmf")
file.parameters["rewrite_function_mesh"] = False
file.parameters["functions_share_mesh"]  = True     # Allows saving of multiple functions
file.parameters["flush_output"]          = True     # Saves the file in case of interruption

# Looping with the update of body force
for (i, inc) in enumerate(Body):
    print("Increment of Body Force: ", Body[i])

    # Solve the problem for the specific step
    solver.solve()

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
