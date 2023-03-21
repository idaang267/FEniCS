# Contact Mechanics Problem using PETSc's TAO solver for nonlinear
# (bound-constrained) optimisation problems

# This example considers a heavy hyperelastic circle under body force in a box
# of the same size

# Copyright (C) 2012 Corrado Maurini
# Modified by Tianyi Li 2014
#
# First added:  2012-09-03
# Last changed: 2014-07-19

# Import Modules
from dolfin import *
import matplotlib.pyplot as plt

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
solver.parameters["report"] = True

# Uncomment this line to see the available parameters
#info(parameters, True)

# Parse (PETSc) parameters
parameters.parse()

# Read pre-made mesh
mesh = Mesh("circle_xyplane.xml")
# Create function space
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

# Kinematics
d = len(u)                       # Spatial dimension
I = Identity(d)                  # Identity tensor
F = I + grad(u)                  # Deformation gradient
C = F.T*F                        # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
B  = Constant((0.0, -1.5))       # Body force per unit volume
E = 10.0                                        # Young's Modulus
nu = 0.3                                        # Poisson Ratio
mu = Constant(E/(2*(1 + nu)))                   # Lamé Parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))    # Lamé Parameter

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic-2 - 2*ln(J)) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - dot(B,u)*dx

# Compute first variation of Pi (directional derivative about u in the
# direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Symmetry condition (to block rigid body rotations)
def symmetry_line(x, on_boundary):
    return near(x[0], 0)
# Subclass funciton space to constrain x (0) direction
bc = DirichletBC(V.sub(0), Constant(0.0), symmetry_line)

# The displacement u must be such that the current configuration
# doesn't escape the box [xmin, xmax] x [ymin, ymax]
constraint_l = Expression(("xmin-x[0]", "ymin-x[1]"), \
                            xmin=-1.0, ymin=-1.0, degree=1)
    # Identifies the left side and bottom side of the square
constraint_u = Expression(("xmax-x[0]", "ymax-x[1]"), \
                            xmax=1.0, ymax=1.0, degree=1)
    # Identifies the right side and top sides of the square
u_min = interpolate(constraint_l, V)
u_max = interpolate(constraint_u, V)

# BC will be incorporated into the lower (min) and upper (max) bounds
bc.apply(u_min.vector())
bc.apply(u_max.vector())

# Solve the problem
solver.solve(ContactProblem(), u.vector(), u_min.vector(), u_max.vector())

# Save solution in XDMF format if available
soln = XDMFFile(mesh.mpi_comm(), "displacement_tao.xdmf")
if has_hdf5():
    soln.write(u)
elif MPI.size(mesh.mpi_comm()) == 1:
    encoding = XDMFFile.Encoding_ASCII
    soln.write(u, encoding)
else:
    # Save solution in vtk format
    soln = File("u.pvd")
    soln << u

# Plot the current configuration
#plot(u, mode="displacement", wireframe=True, title="Displacement field")
#plt.show()
