# Swelling of Unit Cube
#------------------------------------------------------------------------------
# Based on the formulation in 2015 paper "Effect of solvent diffusion on
# crack-tip fields and driving force for fracture of hydrogels" which simplifies
# the free energy formulation in a prior 2015 paper, "A nonlinear, transient FE
# method for coupled solvent diffusion and large deformation of hydrogels"

from dolfin import *                    # Dolfin module
import matplotlib.pyplot as plt         # Module matplotlib for plotting
from mshr import *
from ufl import cofac, rank

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 30,
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
        values[3] = 0.0
    def value_shape(self):
         return (4,)

class OnBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Model parameters
#------------------------------------------------------------------------------
# Name of file
name = "sphere_surface_tension.xdmf"
mesh_res = 10                   # Mesh resolution
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary
chi = 0.6                       # Flory Parameter
l0 = 1.4                        # Initial Stretch (lambda_o)
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Global stepping, chemical stepping, and surface stepping parameters
steps = 0                       # Steps (updated within loop)
tot_steps = 100                 # Total number of time steps
c_steps = 0                     # Chemical step counter (updated within loop)
t_c_steps = 6                   # Total chemical steps
g_steps = 0                     # Surface parameter counter (updated within loop)
t_g_steps = 50                  # Total surface parameter (gamma) steps

# Time parameters
dt = 10**(-3)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression("dt", dt=dt, degree=0)
# Initial time for paraview file
t = 0.0
c_exp = 1.1                     # Control the time step increase
# Surface Energy: Gamma term
gamma = 0.001
Gamma = Expression("gamma", gamma=gamma, degree=0)

# Define mesh
#------------------------------------------------------------------------------
# Create spherical domain with radius of 1.0
domain = Sphere(Point(0.0, 0.0, 0.0), 1.0)
# Discretize mesh by mesh_res
mesh = generate_mesh(domain, mesh_res)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())

# Create boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Measures/redefinition for dx and ds according to subdomains and boundaries
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Define mixed function space
#------------------------------------------------------------------------------
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
# Define Taylor-Hood Elements
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
    # Second order quadratic interpolation for displacement (u)
P1 = FunctionSpace(mesh, "Lagrange", 1)
    # First order linear interpolation for chemical potential
# Use ufl_element() to return the UFL element of the function spaces
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
# Define mixed function space specifying underying finite element
V = FunctionSpace(mesh, P2elem * P1elem)

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

# Define terms for surface tension
N = FacetNormal(mesh)                    # Normal vector in the reference configuration
NansonOp = (cofac(F))                    # Element of area transformation operator
deformed_N = dot(NansonOp,N)             # Surface element vector in the deformed configuration
Jsurf = sqrt(dot(deformed_N,deformed_N)) # Norm of the surface element vector in the current configuration

Isurf = I - outer(N,N)
Fsurf = dot(F, Isurf)

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for surface
# Note: ParaView default view, x is right, y is up, and z-direction is out-of-plane

# Chemical potential BC ramped from mu0 to 0 in the IC class
chem_max = ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3)
chem_p = Expression(("chem_max - (chem_max*c_steps)/t_c_steps"), \
                    chem_max=chem_max, c_steps=c_steps, t_c_steps=t_c_steps, degree=1)

#u_0 = Expression(("(l0-1)*x[0]"), l0=l0, degree=2)

# The Dirichlet BCs are specified in respective subspaces
# Roller displacement BCs
#bc_0 = DirichletBC(V.sub(0).sub(0), u_0, OnBoundary())

# Chemical potential
bc_chem = DirichletBC(V.sub(1), chem_p, OnBoundary())

# Combined boundary conditions
bc = [bc_chem]

# Definitions
#------------------------------------------------------------------------------
# Normalized nominal stress tensor: P = dU/dF
def P(u, mu):
    return F + (-1/J + (1/n)*(1/J + ln((J-1)/J) + chi/(J**2) - mu))*J*inv(F.T)

def Psurf(u):
    return 0.1*inv(Fsurf.T)*Jsurf

# Normalized flux
def Flux(u, mu):
    p1 = dot(inv(F), grad(mu))
    return -(J-1)*dot(p1,inv(F))

# Variational problem where we have two equations for the weak form
F0 = inner(P(u, mu), grad(v_u))*dx - inner(T, v_u)*ds - dot(B, v_u)*dx
F1 = (1/n)*((J-1)*v_mu*dx - (J0-1)*v_mu*dx - DT*dot(Flux(u, mu), grad(v_mu))*dx)
surface_energy_density = Gamma*Jsurf
surface_energy = surface_energy_density*ds
F2 = derivative(surface_energy, u, v_u)
WF = F0 + F1 + F2

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = derivative(WF, w, du)

# SNES solver > Setup Non-linear variational problem
problem = NonlinearVariationalProblem(WF, w, bc, J=Jacobian)
solver_problem = NonlinearVariationalSolver(problem)
solver_problem.parameters.update(snes_solver_parameters)

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile(name)

# Solve for each value using the previous solution as a starting point
while (steps < tot_steps):
    # Print outs to track code progress
    print("Steps: " + str(steps))
    #print("Time: " + str(t))

    # Update fields containing u and mu and solve using the setup parameters
    w0.vector()[:] = w.vector()
    solver_problem.solve()

    # Update the chemical steps for ramping of chemical potential
    if c_steps < t_c_steps:
        c_steps += 1
    c_steps += 0
    chem_p.c_steps = c_steps        # Update steps in expression class

    if g_steps < t_g_steps:
        g_steps += 1
        gamma += 0.001
    Gamma.gamma = gamma

    steps += 1                      # Update total steps

    # Update time parameters
    dt = dt*c_exp;                # Update time step with exponent value
    DT.dt = dt                    # Update time step for weak forms
    t += dt                       # Update total time for paraview file

    # Deep copy not a shallow copy like split(w) for writing the results to file
    (u, mu) = w.split()
    # Project nominal stress to tensor function space
    PTensor = project(P(u, mu), TT)
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential", "mu")
    PTensor.rename("Nominal Stress", "P")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write multiple time-dependent parameters to .xdmf results file
    file_results.write(u,t)
    file_results.write(mu,t)
    file_results.write(PTensor,t)