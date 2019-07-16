# Contact mechanics problem using SNES solver for variational inequalities.

# This example considers a heavy hyperelastic circle under body force in a box
# of the same size

# Copyright (C) 2012 Corrado Maurini
# Modified by Corrado Maurini 2013

# Import Modules
from dolfin import *
import matplotlib.pyplot as plt

# Check if PETSc is installed otherwise exit
if not has_petsc():
    print("DOLFIN must be compiled with PETSc to run this demo.")
    exit(0)

# Define the solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

# Read pre-made mesh with format .xml
mesh = Mesh("circle_xyplane.xml")
# Create vector function space
V = VectorFunctionSpace(mesh, "Lagrange", 1)

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

# Elasticity parameters
B = Constant((0.0, -1.1))                       # Body force per unit volume
E = 10.0                                        # Young's Modulus
nu = 0.3                                        # Poisson Ratio
mu = Constant(E/(2*(1 + nu)))                   # Lamé Parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))    # Lamé Parameter

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 2 - 2*ln(J)) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - dot(B, u)*dx

# Compute first variation of Pi (directional derivative about u in the
# direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# Symmetry condition (to block rigid body rotations)
def symmetry_line(x):
    return abs(x[0]) < DOLFIN_EPS
# Subclass funciton space to constrain x (0) direction
bc = DirichletBC(V.sub(0), 0.0, symmetry_line, method="pointwise")

# The displacement u must be such that the current configuration x+u
# remains in the box [xmin,xmax] x [ymin,ymax]
constraint_l = Expression(("xmin - x[0]","ymin - x[1]"), \
                            xmin=-1.0-DOLFIN_EPS, ymin=-1.0, degree=1)
    # Identifies the left side and bottom side of the square
constraint_u = Expression(("xmax - x[0]", "ymax - x[1]"), \
                            xmax=1.0+DOLFIN_EPS,  ymax=1.0, degree=1)
    # Identifies the right side and top sides of the square
umin = interpolate(constraint_l, V)
umax = interpolate(constraint_u, V)

# Set up the non-linear problem
problem = NonlinearVariationalProblem(F, u, bc, J=J)
problem.set_bounds(umin, umax)

# Set up the non-linear solver
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)
info(solver.parameters, False)

# Solve the problem
(iter, converged) = solver.solve()
print(iter)
# Check for convergence
if not converged:
    warning("This demo is a complex nonlinear problem. Convergence is not guaranteed when modifying some parameters or using PETSC 3.2.")

# Save solution in XDMF format
file = XDMFFile("displacement_snes.xdmf")
file.write(u, 0.)

# plot the current configuration
#plot(u, mode="displacement", wireframe=True, title="Displacement field")
#plt.show()
