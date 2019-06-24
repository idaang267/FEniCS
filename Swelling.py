# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Model is based on Bouklas 2015 Paper "A nonlinear, transient FE method for
# coupled solvent diffusion and large deformation of hydrogels"

from dolfin import *                # Dolfin module
import matplotlib.pyplot as plt     # Module matplotlib

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 200,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "error_on_nonconvergence": False}}

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4

# Initial condition (IC) class for displacement and chemical potential
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1], values[2])
        values[0] = (l0-1)*x[0]
        values[1] = (l0-1)*x[1]
        values[2] = (l0-1)*x[2]
        # Initial Chemical potential: mu0
        values[3] = (1/l0**3)*(l0/bulkMod-1/bulkMod+1) + ln((l0**3-1)/l0**3) + chi/l0**6
    def value_shape(self):
         return (4,)

class Constraint(SubDomain):
    # Top boundary is "target domain" G without one point at [0, 1, 0]
    def inside(self, x, on_boundary):
        return bool((x[1] > 1.0 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS) and \
                    (x[1] > 1.0 - DOLFIN_EPS and x[2] > 0.5 - DOLFIN_EPS) and \
                    on_boundary)
        #return bool([0.5,1,0] and [1,1,0] and [1,1,0.5] and [0.5,1,0.5] and \
        #            [0,1,0.5] and [0,1,1] and [0.5,1,1] and [1,1,1])
    # Map a coordinate x in top boundary to a coordinate y in target boundary (G)
    def map(self, x, y):
        y[0] = x[0]         # x
        y[1] = x[1] - 1.0   # y - 1

        #if near(x[0], 1):
        #    y[0] = x[0] - 1.
        #    y[1] = x[1] - 1.

tb = Constraint()

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
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
l0 = 1.55                       # Initial Stretch (lambda_o)
# Normalization Parameters
bulkMod = 10**3                 # Bulk Modulus (K/N k_B T)
n = 10**(-3)                    # Normalization (N Omega)
# Time Parameters
t = 0.0                         # Initial time, will be updated within loop
num_steps = 3                   # Number of time steps
dt = 10e-1                      # Change in time between time steps
Time = num_steps*dt             # Comparison to control amount of steps

# Define mesh and mixed function space
#------------------------------------------------------------------------------
mesh = UnitCubeMesh(2, 2, 2)                        # Unit Cube
# Define Taylor-Hood Elements
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement (u)
    # Second order quadratic interpolation for displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Chemical Potential (mu)
    # First order linear interpolation for chemical potential
TH = MixedElement([P2, P1])
V  = FunctionSpace(mesh, TH, constrained_domain=tb)

# Just a check
#coordinates = mesh.coordinates()[6]
#print(coordinates)

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for each surface of the cube
# Note: these subdomains are set according to ParaView default view where
# x is right, y is up, and z-direction is out-of-plane
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
back =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)

# Bottom of the cube is attached to substrate and has fixed displacement from
# the reference relative to the current configuration
u_bot = Expression(("(scale-1)*x[0]", "0.0*x[1]", "(scale-1)*x[2]"), scale=l0, degree=1)
    # Account for initial configuration: displacement in x & z direction not y
# Chemical potential BC equivalent to mu0 in the IC class
chem_pot = 0#(1/l0**3)*(l0/bulkMod-1/bulkMod+1) + ln((l0**3-1)/l0**3) + chi/l0**6

# The Dirichlet BCs are specified in respective subspaces
# Displacement in first subspace V.sub(0)
bc_u = DirichletBC(V.sub(0), u_bot, bottom)
# Chemical potential in second subspace V.sub(1)
bc_t = DirichletBC(V.sub(1), chem_pot, top)
bc_l = DirichletBC(V.sub(1), chem_pot, left)
bc_r = DirichletBC(V.sub(1), chem_pot, right)
bc_b = DirichletBC(V.sub(1), chem_pot, back)
bc_f = DirichletBC(V.sub(1), chem_pot, front)
bc = [bc_u, bc_t, bc_l, bc_r, bc_b, bc_f]  # All boundary conditions


# Define functions in mixed function space V
#------------------------------------------------------------------------------
du = TrialFunction(V)                       # Incremental trial function
v = TestFunction(V)                         # Test Function
w = Function(V)                             # Current solution for u and mu
w0 = Function(V)                            # Previous solution for u and mu

# Split test functions and unknowns (produces a shallow copy not a deep copy)
(v_u, v_mu) = split(v)
(u, mu) = split(w)                          # Split current
(u0, mu0) = split(w0)                       # Split previous
    # Deep copy syntax: w.split(deepcopy=True)

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Initial conditions are created by using the class defined and then
# interpolating into a finite element space:
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
Ic = tr(CG)
J = det(F)                      # Current time step
J0 = det(F0)

# Definitions
#------------------------------------------------------------------------------
# Normalized nominal stress tensor: P = dU/dF
def P(u, mu):
    return F + (-1/J + (1/n)*(1/J + ln((J-1)/J) + chi/(J**2) - mu))*J*inv(F.T)

# Definition for normalized flux (Don't use J)
def Flux(u, mu):
    p1 = dot(inv(F), grad(mu))
    return -(J-1)*dot(inv(F), p1)

# Variational problem where we have two equations for the weak form
F0 = inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx
F1 = (1/n)*(inner(J-1, v_mu)*dx - inner(J0-1, v_mu)*dx - inner(Flux(u, mu), grad(v_mu))*dt*dx)
WF = F0 + F1

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = NonlinearVariationalProblem(WF, w, bc, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile("results.xdmf")

# Solve for each value using the previous solution as a starting point
while (t < Time):

    t += dt                        # Update time
    w0.vector()[:] = w.vector()    # Update fields containing u and mu

    # Solve
    solver_problem.solve()

    # Note that this is now a deep copy not a shallow copy like split(w)
    (u, mu) = w.split()

    # Write results to file_results defined
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential","mu")

    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True

    # Write to .xdmf results file
    file_results.write(u,t)
    file_results.write(mu,t)
