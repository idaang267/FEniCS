#------------------------------------------------------------------------------
# Swelling of Unit Sphere (Modeled in 3D)
#------------------------------------------------------------------------------
# Based on the formulation in 2015 paper "Effect of solvent diffusion on
# crack-tip fields and driving force for fracture of hydrogels" which simplifies
# the free energy formulation in a prior 2015 paper, "A nonlinear, transient FE
# method for coupled solvent diffusion and large deformation of hydrogels"

from dolfin import *                    # Dolfin module
import matplotlib.pyplot as plt         # Module matplotlib for plotting
import numpy as np
import os
import shutil
import time

from mshr import *
from ufl import cofac, rank

start_time = time.time()

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4

# Solver parameters: Using PETSc SNES solver
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 100,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-9,
                                          "relative_tolerance": 1e-9,
                                          "error_on_nonconvergence": False}}

# Defining Classes
#------------------------------------------------------------------------------
# Full boundary for chemical potential bc
class OnBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Center point for fixed displacement bc
class pinPoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, tol) and near(x[1], 0.0, tol) and near(x[2], 0.0, tol)

# Right boundary point for fixed displacement bcs (x-axis)
class pinPointRight(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0, tol) and near(x[1], 0.0, tol) and near(x[2], 0.0, tol)

# Front and back boundary point for fixed displacement bcs (z-axis)
class pinPointTop(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, tol) and near(x[1], 1.0, tol) and near(x[2], 0.0, tol)

pinPoint = pinPoint()
pinPointRight = pinPointRight()
pinPointTop = pinPointTop()

# Initial condition (IC) class for displacement and chemical potential
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1], values[2])
        values[0] = (l0-1)*x[0]
        values[1] = (l0-1)*x[1]
        values[2] = (l0-1)*x[2]
        # Initial Chemical potential: mu0
        values[3] = chem_ini #ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3) + gamma/l0
    def value_shape(self):
         return (4,)

# Model parameters
#------------------------------------------------------------------------------
# Set the user parameters, can be parsed from command line
parameters.parse()
userpar = Parameters("user")
userpar.add("chi", 0.2)
userpar.add("gamma", 1.0)
userpar.add("l0", 2.5687)
userpar.add("t_g_steps", 0)
userpar.add("eq_steps1", 0)
userpar.add("t_c_steps", 75)
userpar.add("eq_steps2", 0)
userpar.parse()

# Parse from command line if given
tol = 1E-12
gamma = userpar["gamma"]        # Surface Energy: Gamma term
Gamma = Expression("gamma", gamma=gamma, degree=0)
chi = userpar["chi"]            # Flory Parameter
l0 = userpar["l0"]              # Initial Stretch (lambda_o)
n = 10**(-3)                    # Normalization Parameter (N Omega)
# Global stepping, chemical stepping, and surface stepping parameters
steps = 0                       # Steps (updated within loop)
g_steps = 0                     # Surface parameter counter (updated within loop)
c_steps = 0                     # Chemical step counter (updated within loop)
t_g_steps = userpar["t_g_steps"]# Total surface parameter (gamma) steps
t_c_steps = userpar["t_c_steps"]# Total chemical steps
# Number of steps to reach equilibrium for stress or chemical ramping case
eq_steps1 = userpar["eq_steps1"]
eq_steps2 = userpar["eq_steps2"]
# Total number of time steps
tot_steps = t_g_steps + eq_steps1 + t_c_steps + eq_steps2
# Body force per unit volume and Traction force on the boundary
B = Constant((0.0, 0.0, 0.0))
T = Constant((0.0, 0.0, 0.0))
# Initial and maximum value of chemical potentail
chem_ini = ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3) + n*gamma*2/l0
chem_max = ln((l0**3-1)/l0**3) + 1/l0**3 + chi/l0**6 + n*(1/l0-1/l0**3) + n*gamma*2/l0
# Time parameters
dt = 10**(-2)                   # Starting time step
# Expression for time step for updating in loop
DT = Expression("dt", dt=dt, degree=0)
t = 0.0                         # Initial time (updated in loop)
c_exp = 1.5                     # Control the time step increase

# Name of file
name = "3D_quad"
sim_param1 = "_chi_%.1f" % (chi)
sim_param2 = "_gamma_%.2f" % (gamma)
sim_param3 = "_l0_%.2f" % (l0)
sim_param4 = "_steps_%.0f" % (tot_steps)
savedir = name + sim_param1 + sim_param3 + "/"

# Define mesh
#------------------------------------------------------------------------------
# Create spherical domain with radius of 1.0
mesh = Mesh("sphere_double_quad.xml")

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())

# Create boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Create points
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 3)

# Print out mesh for visualization in Paraview
file_mesh = XDMFFile(savedir + "/" + "mesh.xdmf")
file_mesh.write(mesh)

# # Show points of interest
# points.set_all(0)
# pinPoint.mark(points, 1)
# pinPointTop.mark(points, 1)
# pinPointRight.mark(points, 1)
# file_results = XDMFFile(savedir + "/" + "points.xdmf")
# file_results.write(points)

# Measures/redefinition for dx and ds according to subdomains and boundaries
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Define mixed function space
#------------------------------------------------------------------------------
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
KK = FunctionSpace(mesh,'DG',0)
# Define Taylor-Hood Elements
# Second order quadratic interpolation for displacement (u)
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
# First order linear interpolation for chemical potential
P1 = FunctionSpace(mesh, "Lagrange", 1)
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

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Mark Boundary Subdomains for surface
# Note: ParaView default view, x is right, y is up, and z-direction is out-of-plane

# Displacement BC: pinned center to prevent translation
u_0 = Expression(("0.0","0.0","0.0"), degree=0)
u_dir = Expression(("0.0"), degree=0)

# Chemical potential BC ramped from initial to maximum value
chem_p =  Expression(("c_steps*(chem_max-chem_ini)/t_c_steps + chem_ini"), \
                   chem_ini=chem_ini, chem_max=chem_max, c_steps=c_steps, t_c_steps=t_c_steps, degree=1)

# The Dirichlet BCs are specified in respective subspaces
# Fixed center point to prevent translations
bc_0 = DirichletBC(V.sub(0), u_0, pinPoint, method='pointwise')
# Fixed right point and front point to prevent rotations
bc_pin_RY = DirichletBC(V.sub(0).sub(1), u_dir, pinPointRight, method='pointwise')
bc_pin_RZ = DirichletBC(V.sub(0).sub(2), u_dir, pinPointRight, method='pointwise')
bc_pin_FX = DirichletBC(V.sub(0).sub(0), u_dir, pinPointTop, method='pointwise')
bc_pin_FZ = DirichletBC(V.sub(0).sub(2), u_dir, pinPointTop, method='pointwise')
bc_chem = DirichletBC(V.sub(1), chem_p, OnBoundary())

# Combined boundary conditions
bc = [bc_0, bc_pin_RY, bc_pin_RZ, bc_pin_FX, bc_pin_FZ, bc_chem]

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
Isurf = I - outer(N,N)                   # Surface identity
Fsurf = dot(F, Isurf)                    # Surface deformation gradient

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

# Initialize data arrays
# Saves to show the progression of gamma ramping and chemical ramping
data_steps = np.zeros((tot_steps, 4))

# Save results to an .xdmf file since we have multiple fields (time-dependence)
file_results = XDMFFile(savedir + "/" + name + sim_param1 + sim_param3 + sim_param4 + ".xdmf")

# Loop: solve for each value using the previous solution as a starting point
#------------------------------------------------------------------------------
while (steps < tot_steps):
    if MPI.rank(MPI.comm_world) == 0:
        # Print outs to track code progress
        print("Steps: " + str(steps))
        print("Time: " + str(t))

    # Update fields containing u and mu and solve using the setup parameters
    w0.vector()[:] = w.vector()
    solver_problem.solve()

    # Update the surface stress
    if g_steps < t_g_steps:
        g_steps += 1
        gamma += 0.01
    gamma += 0
    Gamma.gamma = gamma

    # Update the chemical potential
    crit1 = t_g_steps + eq_steps1
    crit2 = t_g_steps + eq_steps1 + t_c_steps
    if (steps > crit1 and steps <= crit2) and c_steps < t_c_steps:
        # Reset time step to original time step
        if c_steps == 0:
            dt = 10**(-3)
        c_steps += 1
    c_steps += 0
    chem_p.c_steps = c_steps        # Update steps in expression class

    chem_val = c_steps*(chem_max-chem_ini)/t_c_steps + chem_ini

    # Save data to plot
    data_steps[steps] = np.array([t, steps, gamma, chem_val])

    # Save data
    #---------------------------------------------------------------------------
    # Deep copy not a shallow copy like split(w) for writing the results to file
    (u, mu) = w.split()
    # Project nominal stress to tensor function space
    PTensor = project(P(u, mu), TT)
    CScalar = project((J-1), KK)
    # Rename results for visualization in Paraview
    u.rename("Displacement", "u")
    mu.rename("Chemical Potential", "mu")
    PTensor.rename("Nominal Stress", "P")
    CScalar.rename("Concentration", "CScalar")
    # Parameters will share the same mesh
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    # Write multiple time-dependent parameters to .xdmf results file
    file_results.write(u,t)
    file_results.write(mu,t)
    file_results.write(PTensor,t)
    file_results.write(CScalar,t)

    # Update total steps and time parameters
    steps += 1                    # Update total steps
    dt = dt*c_exp;                # Update time step with exponent value
    DT.dt = dt                    # Update time step for weak forms
    t += dt                       # Update total time for paraview file

    # if MPI.rank(MPI.comm_world) == 0:
    #     # Save some global quantities as a function of the time
    #     np.savetxt(savedir + '/iterations.txt', data_steps)

# Figures
#---------------------------------------------------------------------------
if MPI.rank(MPI.comm_world) == 0:
    fig, axs = plt.subplots(2, sharex = True)
    fig.suptitle("Ramping Sequence")
    axs[0].plot(data_steps[:, 0], data_steps[:, 2], 'k*-')
    axs[1].plot(data_steps[:, 0], data_steps[:, 3], 'k*-')
    plt.xlabel("Time")
    axs[0].set_title('Gamma Ramping')
    axs[1].set_title('Chemical Potential Ramping')
    plt.savefig(savedir + 'Ramping_Time' + sim_param1 + sim_param2 + sim_param3 + sim_param4 + '.pdf' , transparent=True)
    plt.close()

    fig, axs = plt.subplots(2, sharex = True)
    fig.suptitle("Ramping Sequence")
    axs[0].plot(data_steps[:, 1], data_steps[:, 2], 'k*-')
    axs[1].plot(data_steps[:, 1], data_steps[:, 3], 'k*-')
    plt.xlabel("Steps")
    axs[0].set_title('Gamma Ramping')
    axs[1].set_title('Chemical Potential Ramping')
    plt.savefig(savedir + 'Ramping_Steps' + sim_param1 + sim_param2 + sim_param3 + sim_param4 + '.pdf' , transparent=True)
    plt.close()

    plt.figure(1)
    plt.plot(data_steps[:, 1], data_steps[:, 0], 'k-')
    plt.xlabel("Steps")
    plt.ylabel("Time")
    plt.savefig(savedir + 'TimeProgression' + sim_param1 + sim_param2 + sim_param3 + sim_param4 + '.pdf', transparent=True)
    plt.close()

end_time = time.time()
if MPI.rank(MPI.comm_world) == 0:
    print("Total RunTime: " + str(end_time-start_time))
