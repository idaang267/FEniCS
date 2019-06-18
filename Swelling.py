# Swelling of Unit Cube
# Adapted from Hyperelasticity Demo:

# Import modules
from dolfin import *
import matplotlib.pyplot as plt     # For visualization
import numpy as np
import os

# Class representing the intial conditions for displacement and chemical potential
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1], values[2])
        values[0] = (l0-1)*x[0]
        values[1] = (l0-1)*x[1]
        values[2] = (l0-1)*x[2]
        # Chemical potential mu0
        values[3] = (1/l0**3)*(l0/bulkMod-1/bulkMod+1) + ln((l0**3-1)/l0**3) + chi/l0**6
    def value_shape(self):
         return (4,)

# Class for interfacing with the Newton solver
class Equation(NonlinearProblem):
    # Bilinear and linear forms are used to compute the residual and Jacobian
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)
        self.L = L              # Linear form
        self.a = a              # Bilinear form
        self.bcs = bcs          # Boundary conditions
    def F(self, b, x):          # Computes the residual vector "b"
        assemble(self.L, tensor=b, bcs=self.bcs)
    def J(self, A, x):          # Computes the Jacobian matrix "A"
        assemble(self.a, tensor=A, bcs=self.bcs)

# Model parameters
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.4                       # Flory Parameter
l0 = 1.1                       # Initial Stretch (lambda_o)
# Normalization Parameters
bulkMod = 10^3                 # Bulk Modulus (K/N k_B T)
n = 10^(-3)                     # Normalization (N Omega)
# Time Parameters
t = 0.0
num_steps = 5
dt = 10e-5                      # time step

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True

# Define mesh and mixed function space
mesh = UnitCubeMesh(2, 2, 2)                        # Unit Cube
# Define Taylor-Hood Elements
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement (u)
    # Second order quadratic interpolation for displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Chemical Potential (mu)
    # First order linear interpolation for chemical potential
TH = MixedElement([P2,P1])
V  = FunctionSpace(mesh, TH)

# Boundary Conditions
#------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
# Sides of cube
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
back =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)

# Bottom of the cube is attached to substrate and has fixed displacement
u_bot = Expression(("(scale-1)*x[0]", "(scale-1)*x[1]", "0.0"), scale=l0, degree=1)
# No flux on the bottom surface and four sides
chem_pot = (1/l0**3)*(l0/bulkMod-1/bulkMod+1) + ln((l0**3-1)/l0**3) + chi/l0**6
# The Dirichlet BCs are specified in a subspace
bc_u = DirichletBC(V.sub(0), u_bot, bottom)    # Set displacement on bottom surface

# Set chemical potential boundary conditions to be equivalent to mu0
bc_l = DirichletBC(V.sub(1), chem_pot, left)
bc_r = DirichletBC(V.sub(1), chem_pot, right)
bc_b = DirichletBC(V.sub(1), chem_pot, back)
bc_f = DirichletBC(V.sub(1), chem_pot, front)
bc = [bc_u, bc_l, bc_r, bc_b, bc_f]

# Define trial and test functions of space V
#------------------------------------------------------------------
du = TrialFunction(V)           # Incremental trial function
v = TestFunction(V)             # Test Function
w = Function(V)                 # Current solution
w0 = Function(V)                # Previous solution
# Split (produces a shallow copy not a deep copy)
(v_u, v_mu) = split(v)
(u, mu) = split(w)              #w.split(deepcopy=True)
(u0, mu0) = split(w)            #w0.split(deepcopy=True)

# Initial Conditions
#--------------------------------------------------------------------
# Initial conditions are created by using the classes defined and then
# interpolating the initial conditions into a finite element space:
C_0 = Constant(1 - l0**3)
C0 = interpolate(C_0, V.sub(1).collapse())

# Create intial conditions and interpolate
init = InitialConditions(degree=1)
w0.interpolate(init)
w.interpolate(init)

# Kinematics
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
CG = F.T*F                      # Right Cauchy-Green (CG) tensor

# Invariants of deformation tensors
Ic = tr(CG)
J  = det(F)

# Definition for normalized nominal stress tensor: P = dU/dF
def P(u):
    return F + (-1/J + (1/n)*(1/J + ln((J-1)/J) + chi/J**2 - mu))*J*inv(F)

# Definition for normalized flux (Don't use J)
def Flux(u, mu):
    p1 = dot(inv(F), grad(mu))
    return -(1-J)*dot(inv(F),p1) # where C = - dU/dmu = 1 - J

# Variational problem where we have two equations for the weak form
F0 = inner(P(u), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx
F1 = dot(1-J, v_mu)*dx - inner(C0, v_mu)*dx - inner(Flux(u, mu), grad(v_mu))*dt*dx
WF = F0 + F1

# Compute directional derivative about w in the direction of du (Jacobian)
J = derivative(WF, w, du)

# Create nonlinear problem and Newton solver
problem = Equation(J, WF, bc)
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile("results.xdmf")

# Solve for each value using the previous solution as a starting point
T = num_steps*dt
while (t < T):

    # Update time
    t += dt
    # Update fields containing u and mu
    w0.vector()[:] = w.vector()

    # Solve
    solver.solve(problem, w.vector())
    #solve(WF == 0, w, bc, J=J)

    # Note that this is now a deep copy not a shallow copy like split(w)
    (u, mu) = w.split()               # Split the test functions

    # Write results to file_results defined
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential","mu")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u,t)
    file_results.write(mu,t)
