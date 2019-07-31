# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Model is based on Bouklas 2015 Paper "A nonlinear, transient FE method for
# coupled solvent diffusion and large deformation of hydrogels"

from dolfin import *                    # Dolfin module
import matplotlib.pyplot as plt         # Module matplotlib for plotting
import numpy as np

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 75,
                                          "line_search": "bt",
                                          "linear_solver": "lu",
                                          "method": "vinewtonssls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "report": True,
                                          "error_on_nonconvergence": False}}

#NOTE Need to use vinewtonrsls or vinewtonssls methods if bounds are set.

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4

# Defining Classes
#------------------------------------------------------------------------------
# Initial condition (IC) class for displacement and chemical potential
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1], values[2])
        values[0] = (l0-1)*x[0]
        values[1] = (l0-1)*x[1]
        values[2] = (l0-1)*x[2]
        # Initial Chemical potential: mu0
        values[3] = ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
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
        assemble(self.L, tensor=b)
        for bc in self.bcs:
            bc.apply(b)
    def J(self, A, x):          # Computes the Jacobian matrix "A"
        assemble(self.a, tensor=A)
        for bc in self.bcs:
            bc.apply(A)

# Model parameters
#------------------------------------------------------------------------------
name = "Exp.xdmf"   # Name of file
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
l0 = 1.4                        # Initial Stretch (lambda_o)
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Global stepping and chemical stepping parameters
steps = 0                       # Steps (updated within loop)
c_steps = 0                     # Chemical step counter (updated within loop)
t_c_steps = 6                   # Total chemical steps
tot_steps = 200                 # Total number of time steps
# Time parameters
dt = 10**(-5)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression("dt", dt=dt, degree=0)
t = 0.0                         # Initial time for paraview file
c_exp = 1.2                     # Controls time step increase (10%)

# Define mesh and mixed function space
#------------------------------------------------------------------------------
mesh = UnitCubeMesh(6, 2, 6)              # Unit Cube
TT = TensorFunctionSpace(mesh,'DG',0)       # Tensor space for stress projection
V0 = FunctionSpace(mesh, "DG", 0)           # Vector space for contact pressure
# Taylor-Hood Elements for displacment (u) and chemical potential (mu)
# 2nd order quadratic interpolation for u
P2 = VectorFunctionSpace(mesh, "Lagrange", 2, dim = 3)
# 1st order linear interpolation for mu
P1 = FunctionSpace(mesh, "Lagrange", 1)
# Use ufl_element() to return the UFL element for each function space
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
# Define mixed function space specifying underlying finite element
V = FunctionSpace(mesh, P2elem * P1elem)

# Define functions in mixed function space V
#------------------------------------------------------------------------------
du = TrialFunction(V)                       # Incremental trial function
v = TestFunction(V)                         # Test Function
w = Function(V)                             # Current solution for u and mu
w0 = Function(V)                            # Previous solution for u and mu
# Split test functions and unknowns (produces a shallow copy not a deep copy)
(v_u, v_mu) = split(v)                      # Split test function
(u, mu) = split(w)                          # Split current solution
(u0, mu0) = split(w0)                       # Split previous solution

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube
# Note: these subdomains are set according to ParaView default view where
# x is right, y is up, and z-direction is out-of-plane
bot = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
# Lateral surfaces: left and right
def symmetry_x(x, on_boundary):
    return near(x[0], 0) or near(x[0],1.0)
# Lateral surface: back and front
def symmetry_z(x, on_boundary):
    return near(x[2], 0) or near(x[2],1.0)

# Bottom of the cube is attached to substrate and has fixed displacement from
# the reference relative to the current configuration
u_bot = Expression(("(l0-1)*x[0]", "0.0*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
    # Account for initial configuration: displacement in x & z direction not y

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(l0-1)*x[0]"), l0=l0, degree=1)
u_bf = Expression(("(l0-1)*x[2]"), l0=l0, degree=1)
# Chemical potential BC ramped from 0 to mu0 in the IC class
chem_max = ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
chem_p = Expression(("chem_max - (chem_max*c_steps)/t_c_steps"), \
                    chem_max=chem_max, c_steps=c_steps, t_c_steps=t_c_steps, degree=1)

# The Dirichlet BCs are specified in respective subspaces
# Displacement in first subspace V.sub(0)
bc_u = DirichletBC(V.sub(0), u_bot, bot)

# Roller displacement BCs on lateral faces
# Sub-subclass function space to constrain normal degrees of freedom
# V.sub(0).sub(0) gives x direction
bc_r_l = DirichletBC(V.sub(0).sub(0), u_lr, symmetry_x)
bc_r_r = DirichletBC(V.sub(0).sub(0), u_lr, symmetry_x)
# V.sub(0).sub(2) gives z direction
bc_r_b = DirichletBC(V.sub(0).sub(2), u_bf, symmetry_z)
bc_r_f = DirichletBC(V.sub(0).sub(2), u_bf, symmetry_z)

# Chemical potential in second subspace V.sub(1) on top surface
bc_t = DirichletBC(V.sub(1), chem_p, top)

#===============================================================================
# Define the indenter (obstacle) with a parabola
R, d = 0.5, 0.02
#constraint_l = Expression(("x[0]","-d+(pow(x[0],2)+pow(x[2],2))/(2*R)","x[2]"), l0=l0, d=d, R=R, degree=2)
constraint_l = Expression(("-d+(pow(x[0],2)+pow(x[2],2))/(2*R)"), d=d, R=R, degree=2)

lower = Function(V)
upper = Function(V)

ninfty2 = Function(P2)
ninfty1 = Function(P1)
pinfty2 = Function(P2)
pinfty1 = Function(P1)

ninfty2.vector()[:] = -np.infty
ninfty1.vector()[:] = -np.infty
pinfty2.vector()[:] = np.infty
pinfty1.vector()[:] = np.infty

fa = FunctionAssigner(V, [P2, P1])
print(fa.num_receiving_functions()) # V
print(fa.num_assigning_functions()) # [P2, P1]

#fa.assign(lower, [interpolate(constraint_l, P2), ninfty1])
fa.assign(lower, [ninfty2, ninfty1])
fa.assign(upper, [pinfty2, pinfty1])
#===============================================================================

# Combined boundary conditions
bc = [bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f, bc_t]

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Initial conditions are created by using the class defined and then
# interpolating into a finite element space
init = InitialConditions(degree=1)          # Expression requires degree def.
w.interpolate(init)                         # Interpolate current solution
w0.interpolate(init)                        # Interpolate previous solution

# Kinematics
#------------------------------------------------------------------------------
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient from current time step
F0 = I + grad(u0)               # Deformation gradient from previous time step
CG = F.T*F                      # Right Cauchy-Green (CG) tensor

# Invariants of deformation tensors
Ic = tr(CG)                     # First invariant
J = det(F)                      # Current time step for third invariant
J0 = det(F0)                    # Previous time step for third invariant

# Definitions
#------------------------------------------------------------------------------
# Normalized nominal stress tensor: P = dU/dF
def P(u, mu):
    return F + (-1/J + (1/n)*(1/J + ln((J-1)/J) + chi/(J**2) - mu))*J*inv(F.T)

# Normalized flux
def Flux(u, mu):
    p1 = dot(inv(F), grad(mu))
    return -(J-1)*dot(p1,inv(F))

# Variational problem where we have two equations for the weak form
F0 = inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx
F1 = (1/n)*((J-1)*v_mu*dx - (J0-1)*v_mu*dx - DT*dot(Flux(u, mu), grad(v_mu))*dx)
WF = F0 + F1

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = NonlinearVariationalProblem(WF, w, bc, J=Jacobian)
problem.set_bounds(lower, upper)

solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile(name)

# Define contact outputs (names for Paraview)
p = Function(V0, name="Contact pressure")

# Solve for each value using the previous solution as a starting point
while (steps < tot_steps):

    # Update the chemical steps for ramping of chemical potential
    if c_steps < t_c_steps:
        c_steps += 1
    c_steps += 0
    chem_p.c_steps = c_steps        # Update steps in expression class

    # Update fields containing u and mu and solve using the setup parameters
    w0.vector()[:] = w.vector()
    solver.solve()

    steps += 1                      # Update total steps

    # Note that this is now a deep copy not a shallow copy like split(w)
    # allowing us to write the results to file
    (u, mu) = w.split()
    # Project nominal stress
    PTensor = project(P(u, mu), TT)

    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential", "mu")
    PTensor.rename("Nominal Stress", "P")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(u,t)
    file_results.write(mu,t)
    file_results.write(PTensor,t)
    file_results.write(p, t)

    # Update time parameters
    dt = dt*c_exp;                # Update time step with exponent value
    DT.dt = dt                    # Update time step for weak forms
    t += dt                       # Update total time for paraview file

    # Print to track code progress
    print(steps)
