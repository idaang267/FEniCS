# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Model is based on Bouklas 2015 Paper "A nonlinear, transient FE method for
# coupled solvent diffusion and large deformation of hydrogels"

from dolfin import *                    # Dolfin module
import matplotlib.pyplot as plt         # Module matplotlib for plotting

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 100,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          # Newton line search
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "error_on_nonconvergence": False}}

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
        values[3] = 0 #ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
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
name = "penalty.xdmf"   # Name of file
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
l0 = 1.4                        # Initial Stretch (lambda_o)
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Global stepping and chemical stepping parameters
steps = 0                       # Steps (updated within loop)
c_steps = 0                     # Chemical step counter (updated within loop)
t_c_steps = 6                   # Total chemical steps
t_indent_steps = 5              # Total indentation steps
tot_steps = 10                  # Total number of time steps
# Time parameters
dt = 10**(-5)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression ("dt", dt=dt, degree=0)
t = 0.0                         # Initial time for paraview file
c_exp = 1.2                     # Controls time step increase (10% )
# Define mesh
N_plane = 16                               # Number of elements on top plane
mesh = UnitCubeMesh(N_plane, 3, N_plane)   # Unit Cube

# Define function spaces
#------------------------------------------------------------------------------
TT = TensorFunctionSpace(mesh,'DG',0)      # Tensor space for stress projection
V0 = FunctionSpace(mesh, "DG", 0)          # Vector space for contact pressure
# Taylor-Hood Elements for displacment (u) and chemical potential (mu)
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
    # Second order quadratic interpolation for u
P1 = FunctionSpace(mesh, "Lagrange", 1)
    # First order linear interpolation for mu
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
# Lateral surfaces
left  = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
back  = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[2], side) && on_boundary", side = 1.0)

# Exterior facets
#------------------------------------------------------------------------------
# Note, see 1.4.0 #24 "Marking subdomains of a mesh" for more information
# Use MeshFunction() to store the numbering of the subdomains. An argument
# specifying the type of the MeshFunction must be given, where allowed types
# are ‘int’, ‘size_t’, ‘double’ and ‘bool’. Second argument (optional) specifies
# the mesh, while the third argument (optional) gives the topological dimension
facets = MeshFunction("size_t", mesh, 2)
# Mark all facets as part of subdomain 0
facets.set_all(0)
# Top exterior facets are marked as subdomain 1, using 'top' boundary
top.mark(facets, 1)
# Measure redefines ds
ds = Measure('ds', subdomain_data=facets)

# Bottom of the cube is attached to substrate and has fixed displacement from
# the reference relative to the current configuration
u_bot = Expression(("(l0-1)*x[0]", "0.0*x[1]", "(l0-1)*x[2]"), l0=l0, degree=1)
    # Account for initial configuration: displacement in x & z direction not y

# Roller boundary conditions (define normal directions) on lateral surfaces
u_lr = Expression(("(l0-1)*x[0]"), l0=l0, degree=1)
u_bf = Expression(("(l0-1)*x[2]"), l0=l0, degree=1)
# Chemical potential BC ramped from 0 to mu0 in the IC class
# Temp set to 0
chem_max = 0#ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
chem_p = Expression(("chem_max - (chem_max*c_steps)/t_c_steps"), \
                    chem_max=chem_max, c_steps=c_steps, t_c_steps=t_c_steps, degree=1)

# The Dirichlet BCs are specified in respective subspaces
# Displacement in first subspace V.sub(0)
bc_u = DirichletBC(V.sub(0), u_bot, bot)

# Roller displacement BCs on lateral faces
# Access normal degree of freedom (Ex. V.sub(0).sub(0) gives x direction)
bc_r_l = DirichletBC(V.sub(0).sub(0), u_lr, left)
bc_r_r = DirichletBC(V.sub(0).sub(0), u_lr, right)
bc_r_b = DirichletBC(V.sub(0).sub(2), u_bf, back)
bc_r_f = DirichletBC(V.sub(0).sub(2), u_bf, front)

# Chemical potential in second subspace V.sub(1) on top surface
bc_t = DirichletBC(V.sub(1), chem_p, top)

# Combined boundary conditions
bc = [bc_u, bc_r_l, bc_r_r, bc_r_b, bc_r_f]#, bc_t]

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

# Definition of The Mackauley bracket <x>+
def ppos(x):
    return (x+abs(x))/2.

# Define the indenter with a parabola
R = 0.25                        # Indenter radius
depth = 0.1                     # Initial indenter depth of indentation
indent = Expression("-d+l0-1+(pow((x[0]-0.5),2)+pow((x[2]-0.5),2))/(2*(R+l0-1))", \
                        l0=l0, d=depth, R=R, degree=2)

# Note: A large penalty parameter deteriorates the problem conditioning,
# solving time will drastically increase and the problem can fail
pen = Constant(1e4)

# Variational problem where we have two equations for the weak form
F0 = inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx
F1 = (1/n)*((J-1)*v_mu*dx - (J0-1)*v_mu*dx - DT*dot(Flux(u, mu), grad(v_mu))*dx)
# Penalty method, access the y[1] direction and call subdomain 1 (top)
F2 = pen*dot(v_u[1], ppos(u[1]-indent))*ds(1)
WF = F0 + F1 + F2

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = NonlinearVariationalProblem(WF, w, bc, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

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
    chem_p.c_steps = c_steps                # Update steps in expression class

    # Update indenter indentation depth
    if t_indent_steps > steps:
        depth += 0.05
        indent.depth = depth                # Update depth in expression class

    # Update fields containing u and mu and solve using the setup parameters
    w0.vector()[:] = w.vector()
    solver_problem.solve()

    steps += 1                              # Update total steps

    # Note that this is now a deep copy not a shallow copy like split(w)
    # allowing us to write the results to file
    (u, mu) = w.split()
    PTensor = project(P(u, mu), TT)         # Project nominal stress
    p.assign(-project(P(u, mu)[1, 1], V0))  # Project pressure

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
    print("Radius: " + str(R))
    print("Depth: " + str(depth))
