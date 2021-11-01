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
from ufl import rank

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"
# The following two fields were in the original hyperelasticity demo
parameters["form_compiler"]["cpp_optimize"] = True
# parameters["form_compiler"]["representation"] = "uflacs"
    # For processing UFLACS - unified form language expressions
# parameters["form_compiler"]["log_level"] = INFO
    # Show progress of compiler
# parameters["allow_extrapolation"] = True
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

# Element-wise projection using LocalSolver
def local_project(v, V, u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

# Sub domain for clamp at left end
def left(x, on_boundary):
    return near(x[0], 0.) and on_boundary

# Sub domain for rotation at right end
def right(x, on_boundary):
    return near(x[0], 1.) and on_boundary

# Control the number of loop iterations through user parameters
parameters.parse()
user_par = Parameters("user")
user_par.add("p0", 1.0)
user_par.add("T", 3.0)
user_par.add("Nsteps", 50)      # Displacement step number
user_par.parse()

# Elasticity parameters
E, nu = 1000.0, 0.3
mu    = Constant(E/(2.0*(1.0 + nu)))                  # Lamé Parameter
lmbda = Constant(E*nu/((1.0 + nu)*(1.0 - 2.0*nu)))    # Lamé Parameter

# Time-discretization parameters: Generalized-alpha method parameters
# alpha_m=0.2` and alpha_f=0.4 ensures unconditional stability
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
# Optimal dissipation and second-order accuracy choice for beta and gamma
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)

# Define the time-stepping parameters
T   = user_par["T"]      # Final time of the interval
Nsteps = user_par["Nsteps"]    # Number of time steps
dt = Constant(T/Nsteps)     # compute associated time interval between two steps
p0 = user_par["p0"]
cutoff_Tc = T/2
# Define the loading of JIT-compiled expression depending on t using conditional syntax
p = Expression(("0", "t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)

# Definitions of constants from hyperelasticity demo
B = Constant((0.0, 0.0, 0.0))  # Body force per unit volume

# Naming parameters
simulation_params = "S_%.0f" % (Nsteps)
savedir   = "output/"+simulation_params+"/"

# Define mesh and mixed function space
mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 60, 10, 5)

# Create mesh function over the cell facets for distinguishing different boundaries
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
# Mark the right extremity
force_boundary = AutoSubDomain(right)
force_boundary.mark(boundary_subdomains, 3)

# Mark boundary subdomains
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

# Define measure for boundary condition integral
# The exterior surface measure, ds, defined using the boundary subdomains
dss = ds(subdomain_data=boundary_subdomains)

# Taylor-Hood space
# P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement
# P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Hydrostatic Pressure
# TH = MixedElement([P2,P1])
# V  = FunctionSpace(mesh, TH)
V = VectorFunctionSpace(mesh, "CG", 1)
# Define tensorial DG-0 function space for saving stress field evolution
Vsig = TensorFunctionSpace(mesh, "DG", 0)

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, hydrostatic pressure)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
# Functions including the solution (u) from previous iteration
u  = Function(V, name="Displacement")
Stress = Function(Vsig, name="sigma")

# Fields from previous time step (displacement, velocity, acceleration) from
# previous increment t_n
u_old = Function(V)
v_old = Function(V)
a_old = Function(V)

# Define Dirichlet boundary (x = 0 or x = 1)
c = Constant((0.0, 0.0, 0.0))
# The Dirichlet BCs are specified in a subspace
bcl = DirichletBC(V, c, left)
bcs = [bcl]

# Kinematics
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Where P = dW/dF:
def P(u):
    return mu*(F - inv(F.T))

def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new

psi = (mu/2)*(Ic - 3 - 2*ln(J))
Pi = psi*dx

# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
a_new = (u-u_old-dt*v_old)/(beta*dt**2) - (1-2*beta)/(2*beta)*a_old
# Update formula for velocity: v = dt * ((1-gamma)*a0 + gamma*a) + v0
v_new = v_old + dt*((1-gamma)*a_old + gamma*a_new)

# Specify the quadrature degree for efficiency
WF = inner(P(avg(u_old, u, alpha_f)), grad(v))*dx - dot(B, v)*dx - dot(p, v)*dss(3) \
   + alpha_f*dot(a_old, v)*dx + (1-alpha_f)*dot(a_new, v)*dx

# Directional derivative
J_o = derivative(WF, u, du)

# Solve variational problem
varproblem = NonlinearVariationalProblem(WF, u, bcs, J=J_o)
solver = NonlinearVariationalSolver(varproblem)
solver.parameters.update(snes_solver_parameters)
# info(solver.parameters, True)
# (iter, converged) = solver.solve()

# Loading parameter (list of values of surface tension for the simulations)
time = np.linspace(0, T, Nsteps+1)      # Time-stepping
u_tip = np.zeros((Nsteps+1,))           # Displacement at tip of sample

# Save results to an .xdmf file since we have multiple fields
# Save results to an .xdmf file since we have multiple fields
xdmf_file = XDMFFile(savedir + "elastodynamics-results.xdmf")
# Parameters will share the same mesh
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = False

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
for (i, dt) in enumerate(np.diff(time)):

    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    t = time[i+1]

    # Structure used for one printout of the statement
    print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i, t))

    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    p.t = t - float(alpha_f*dt)

    # Solve the nonlinear problem (using Newton solver)
    solver.solve()

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()
    # use update functions using vector arguments
    a_vec = (u_vec-u0_vec-dt*v0_vec)/(beta*dt**2) - (1-2*beta)/(2*beta)*a0_vec
    v_vec = v0_vec + dt*((1-gamma)*a0_vec + gamma*a_vec)
    # Update (u_old <- u)
    u_old.vector()[:] = u.vector()
    v_old.vector()[:] = v_vec
    a_old.vector()[:] = a_vec

    local_project(P(u), Vsig, Stress)
    local_project(v_old, V, v_old)
    local_project(a_old, V, a_old)

    v_old.rename("Velocity", "v_old")
    a_old.rename("Acceleration", "a_old")

    xdmf_file.write(u,t)
    xdmf_file.write(Stress, t)
    xdmf_file.write(v_old, t)
    xdmf_file.write(a_old, t)

    p.t = t

    # Record tip displacement and compute energies
    u_tip[i+1] = u(1., 0.05, 0.)[1]

# Plot tip displacement evolution
plt.figure()
plt.plot(time, u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
plt.ylim(-0.5, 0.5)
plt.savefig(savedir + "/DisplacementEvolution.pdf", transparent=True)
plt.close()
