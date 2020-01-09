# =============================================================================
# FEniCS code  Variational Fracture Mechanics
# =============================================================================
#
# A static solution of the variational fracture mechanics problems using the
# regularization two-fold anisotropic damage model
#
# author: Bin Li (bl736@cornell.edu), Corrado Maurini (corrado.maurini@upmc.fr)
# Edited: Ida Ang & Jason Mulderrig
#
# date: 10/10/2017, 12/17/2019

# =============================================================================
# Import all packages from dolfin and mshr
from dolfin import *
from mshr import *

# Import specific packages
from ufl import replace
import math
import os
import shutil
import sympy
import numpy as np
import matplotlib.pyplot as plt

# Define the minimization problem by using OptimisationProblem class
# ----------------------------------------------------------------------------
class DamageProblem(OptimisationProblem):
    # Constructor
    def __init__(self,f,gradf,J,alpha,bcs):
        OptimisationProblem.__init__(self)
        self.total_energy = f
        self.Dalpha_total_energy = gradf
        self.J_alpha = J
        self.alpha = alpha
        self.bc_alpha = bcs
    # Weak Form
    def f(self, x):
        self.alpha.vector()[:] = x
        return assemble(self.total_energy)
    # First directional derivative: First variation
    def F(self, b, x):
        self.alpha.vector()[:] = x
        assemble(self.Dalpha_total_energy, b)
        for bc in self.bc_alpha:        # Apply boundary conditions
            bc.apply(b)
    # Gateaux derivative: Second variation
    def J(self, A, x):
        self.alpha.vector()[:] = x
        assemble(self.J_alpha, A)
        for bc in self.bc_alpha:        # Apply boundary conditions
            bc.apply(A)

# Parameters for DOLFIN and SOLVER
# ----------------------------------------------------------------------------
# log level: LogLevel.WARNING or INFO
set_log_level(LogLevel.WARNING)
# set some dolfin specific parameters
info(parameters,True)
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

# Parameters of the linear solver for displacement (u) problem
# -----------------------------------------------------------------------------
solver_u_parameters = {"linear_solver": "mumps",
                       "symmetric": True,
                       "preconditioner": "hypre_amg",
                       "krylov_solver": {"report": False,
                                         "monitor_convergence": False,
                                         "relative_tolerance": 1e-10,
                                         "absolute_tolerance": 1e-10}}

# Parameters of the PETSc/Tao solver used for the damage (alpha) problem
tao_solver_parameters = {"maximum_iterations": 200,
                         "report": False,
                         "line_search": "more-thuente",
                         "linear_solver": "mumps",
                         "method": "tron",
                         "gradient_absolute_tol": 1e-8,
                         "gradient_relative_tol": 1e-8,
                         "error_on_nonconvergence": False}

# Set the user parameters
# ------------------------------------------------------------------------------
# Can specify parameters using python3 _.py --meshsize 300
parameters.parse()
userpar = Parameters("user")
userpar.add("meshsize",200)
userpar.add("load_min",0.)
userpar.add("load_max",4.0)
userpar.add("load_steps",201)
userpar.parse()

# Parameters for ANISOTROPIC surface energy and materials
# ------------------------------------------------------------------------------
# Material constants
E     = Constant(1.0)   # Young's Modulus
nu    = Constant(0.3)   # Poisson's Ratio
Gc    = Constant(1.0)   # Griffith's criterion: energy release rate
k_ell = Constant(1.e-6) # Residual stiffness

# Loading Parameters
ut    = 1.0             # Reference value for the loading (imposed displacement)
w1    = 1.0
E_0   = 1.0

# Numerical parameters of the alternate minimization
maxiteration = 2000
AM_tolerance = 1e-4

# Geometry paramaters
L = 1.0
H = 0.1
N = userpar["meshsize"]
# Used as a tolerance value to find boundaries
hsize = float(L/N)
# damage parameter: internal length scale used for tuning Gc
ell   = Constant(6.0*hsize)

# Naming parameters for saving output
case = "B"
modelname = "traction-bar"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "h_%.4f" % (hsize)
savedir   = "output/" + case + "_" + simulation_params + "/"
filename = "Example" + case

# If directory exists, remove recursively and create new directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Last parameter of RectangleMesh is the direction of diagonals: alternating diagonals
mesh = RectangleMesh(Point(0., -0.5*H), Point(L, 0.5*H), int(N), int(float(H/hsize)), "right/left")
# Save mesh as XDMF file for visualization in Paraview
geo_mesh = XDMFFile(MPI.comm_world, savedir+meshname)
geo_mesh.write(mesh)

# First call to mesh only creates entities of dimension zero (vertices) and
# entities of the maximal dimension (cells). Other entities must be explicitly
# created by calling init()
mesh.init()
# Obtain number of space dimensions (1D, 2D, 3D)
ndim = mesh.geometry().dim()
# Print dimensions of mesh
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# Strain, Stress, and Constitutive functions of the damage model
# ----------------------------------------------------------------------------
def eps(v):                   # Strain
    return sym(grad(v))
def sigma_0(v):               # Stress following Hooke's Law (Plane stress)
    mu = E/(2.0*(1.0+nu))
    lmbda = E*nu/(1.0-nu**2)
    return 2.0*mu*(eps(v))+lmbda*tr(eps(v))*Identity(ndim)

# Constitutive functions of the damage model
def w(alpha):
    w1 = 1.0
    # Depending on the model type, alpha or alpha**2
    if case == "A":
        return w1*alpha
    else:
        return w1*alpha**2

def a(alpha):
    return (1.0-alpha)**2

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
def left_boundary(x, on_boundary):
    return on_boundary and near(x[0], 0.0, 0.1 * hsize)
def right_boundary(x, on_boundary):
    return on_boundary and near(x[0], L, 0.1 * hsize)

# Note when using "pointwise", the boolean argument on_boundary in SubDomain
# will always be false
def left_pinpoints(x, on_boundary):     # center left point
    return near(x[0], 0.0, 0.1 * hsize) and near(x[1], 0.0, 0.1 * hsize)
def right_pinpoints(x, on_boundary):    # center right point
    return near(x[0], L, 0.1 * hsize) and near(x[1], 0.0, 0.1 * hsize)

# Variational formulation
# ----------------------------------------------------------------------------
# Create function spaces for 2D elasticity problem and damage problem
V_u     = VectorFunctionSpace(mesh, "Lagrange", 1)
V_alpha = FunctionSpace(mesh, "Lagrange", 1)
# Tensor Function space for stress and strain
TT = TensorFunctionSpace(mesh,'Lagrange', 1)

# Define the function, test and trial fields for elasticity and damage problem
u  = Function(V_u, name="Displacement")
du = TrialFunction(V_u)
v  = TestFunction(V_u)
alpha  = Function(V_alpha, name="Damage")
dalpha = TrialFunction(V_alpha)
beta   = TestFunction(V_alpha)

# For postprocessing, where we want values along a certain section of the mesh
# ----------------------------------------------------------------------------
# Returns number of nodal points for specific function space
dim = V_alpha.dim()
dimTT = TT.dim()
# Degrees of Freedom (dof)
dof = V_alpha.tabulate_dof_coordinates().reshape(dim, ndim)
dofTT = TT.tabulate_dof_coordinates().reshape(dimTT, ndim)
# Cound find coordinates spanning x by dof[:, 0] and z (if 3D) by dof[:, 2]
# Coordinates spanning y
y = dof[:, 1]
yTT = dofTT[:, 1]
# Extract dof indices where some condition is met: here the centerline
indices = np.where(y == 0.0)[0]
indicesTT = np.where(yTT == 0.0)[0]

# Dirichlet (displacement, damage) boundary conditions (BCs)
# --------------------------------------------------------------------
# Impose the displacements field given by asymptotic expansion of crack tip
u_UL = Expression("0.0", degree=0)
# Sliding/Pulling displacement BC: Initialized to 0 and updated in the loop
u_UR = Expression("t", t=0.0, degree=0)

# Boundary conditions - u (displacement)
# Roller BCs: No displacement in the x-direction of the left and right boundary
Gamma_u_0 = DirichletBC(V_u.sub(0), u_UL, left_boundary)
Gamma_u_1 = DirichletBC(V_u.sub(0), u_UR, right_boundary)
# Slider boundary condition
Gamma_u_2 = DirichletBC(V_u.sub(1), u_UL, left_pinpoints, method='pointwise')
#Gamma_u_3 = DirichletBC(V_u.sub(1), u_UL, right_pinpoints, method='pointwise')

# Combine displacement boundary conditions
bc_u = [Gamma_u_0, Gamma_u_1, Gamma_u_2]

# Boundary conditions - alpha (zero damage)
# Apply to space V_alpha. No damage on the left and right boundary
Gamma_alpha_0 = DirichletBC(V_alpha, 0.0, left_boundary)
Gamma_alpha_1 = DirichletBC(V_alpha, 0.0, right_boundary)

# Combine damage boundary conditions
bc_alpha = [Gamma_alpha_0, Gamma_alpha_1]

# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Fenics forms for the energies
def sigma(u, alpha):
    return (a(alpha)+k_ell)*sigma_0(u)

# Zero body force defined
body_force        = Constant((0., 0.))
# Elastic Energy 0.5(Stress : Strain)
elastic_energy    = 0.5*inner(sigma(u, alpha), eps(u))*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Keep for reference because this is familiar to HE demo. We already know that
# this is a linear problem so we don't need to use a nonlinear solver.
'''
# Weak form of elasticity problem
E_u  = derivative(elastic_potential, u, v)
# Variational problem for the displacement
problem_u = NonlinearVariationalSolver(problem_u)
# Set up the solvers
solver_u = NonlinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
'''

# Weak form of elasticity problem
E_u  = derivative(elastic_potential, u, v)
# Writing tangent problems in term of test & trial functions for matrix assembly
E_du = replace(E_u, {u: du}) # Replace u by du in E_u equation
# Variational problem for the displacement
problem_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), u, bc_u)
# Set up the solvers
solver_u  = LinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
# Output to see solver parameters
# info(solver_u.parameters, True)

# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Initialization of field (known) alpha
alpha_0 = interpolate(Expression("0.0", degree=0), V_alpha)
# Declare z as a symbol to do symbolic integration
z = sympy.Symbol("z", positive=True)
# Damage parameter: Normalization constant c_w
# Integrate square root of w(z) with respect to z from bounds 0 to 1
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*inner(grad(alpha), grad(alpha)))*dx
damage_functional = elastic_potential + dissipated_energy

# Compute directional derivative about alpha in the direction of beta (Gradient)
E_alpha = derivative(damage_functional, alpha, beta)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
E_alpha_alpha = derivative(E_alpha, alpha, dalpha)

# Implement box constraints (lower and upper bounds) for damage field
# --------------------------------------------------------------------
# Variational problem for the damage
# Use non-linear variational inequality solvers of PETSc
damage_problem = DamageProblem(damage_functional, E_alpha, E_alpha_alpha, alpha, bc_alpha)

# Set up the PETSc TAO solver
solver_alpha  = PETScTAOSolver()
# Update parameters
solver_alpha.parameters.update(tao_solver_parameters)
# Set an upper and lower bound for damage variable between [0, 1]
alpha_lb = interpolate(Expression("0.0", degree=0), V_alpha)  # lower bound, set to 0
alpha_ub = interpolate(Expression("1.0", degree=0), V_alpha)  # upper bound, set to 1

# Loading and initialization of vectors to store data
# ----------------------------------------------------------------------------
load_multipliers = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies         = np.zeros((len(load_multipliers), 4))
iterations       = np.zeros((len(load_multipliers), 2))

# Stores dimensionalized and nondimensionalized date for each time step
# Order: Strain, damage variable alpha, stress
dim_data    = np.zeros((len(load_multipliers), 3))
nondim_data = np.zeros((len(load_multipliers), 3))

# Set the saved data file name
file_results     = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
file_results.parameters["rewrite_function_mesh"] = False
file_results.parameters["functions_share_mesh"]  = True
file_results.parameters["flush_output"]          = True
# Write the parameters to file
File(savedir+"/parameters.xml") << userpar

# Constants for Examples A/1 and B/2 for Gradient Damage Models
if case == "A":
    sigma_e = sqrt(w1*E_0)
    U_e = (sigma_e/E_0)*L
else:
    sigma_m = (3*sqrt(3))/(8*sqrt(2))*sqrt(w1*E_0)
    U_m = (16/9)*(sigma_m/E_0)*L

# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    # Update time for displacement boundary condition
    u_UR.t = t * ut

    # Print updated time
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))

    # Alternate Mininimization Scheme
    # -------------------------------------------------------------------------
    # Solve for u holding alpha constant then solve for alpha holding u constant
    iteration = 1       # Initialization of iteration loop
    err_alpha = 1.0     # Initialization for condition for iteration
    # Conditions for iteration
    while err_alpha > AM_tolerance and iteration < maxiteration:
        # Solve elastic problem
        solver_u.solve()
        # Solve damage problem with box constraint on alpha to [1, 1]
        solver_alpha.solve(damage_problem, alpha.vector(), alpha_lb.vector(), alpha_ub.vector())
        # Update the alpha condition for iteration by calculating the alpha error norm
        alpha_error = alpha.vector() - alpha_0.vector()
        err_alpha = alpha_error.norm('linf')        # Row-wise norm
        # Monitor the number of iterations
        if MPI.rank(MPI.comm_world) == 0:
            print ("AM Iteration: {0:3d},  alpha_error: {1:>14.8f}".format(iteration, err_alpha))

        # Update iteration
        alpha_0.assign(alpha)
        iteration = iteration + 1

    # Updating the lower bound to account for the irreversibility of damage
    alpha_lb.vector()[:] = alpha.vector()

    # --------------------------------------------------------------------------
    # Project strain and stress to tensor function space
    strain_val = project(eps(u), TT)
    sigma_val = project(sigma(u, alpha), TT)

    # Rename variables for visualization
    alpha.rename("Damage", "alpha")
    u.rename("Displacement", "u")
    strain_val.rename("Strain", "strain_val")
    sigma_val.rename("Stress", "sigma_val")

    # Post-processing
    # -------------------------------------------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value = assemble(elastic_energy)
    surface_energy_value = assemble(dissipated_energy)
    # save energies
    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, elastic_energy_value+surface_energy_value])

    # Save solutions to file for visualization to Paraview
    file_results.write(strain_val, t)
    file_results.write(sigma_val, t)
    file_results.write(alpha, t)
    file_results.write(u, t)

    if MPI.rank(MPI.comm_world) == 0:
        print("\nEnd of timestep {0:3d} with load multiplier {1:4f}".format(i_t, t))
        print("\nElastic and Surface Energies: [{0:6f},{1:6f}]".format(elastic_energy_value, surface_energy_value))
        print("\nElastic and Surface Energies: [{},{}]".format(elastic_energy_value, surface_energy_value))
        print("-----------------------------------------")
        # Save some global quantities as a function of the time
        np.savetxt(savedir + '/energies.txt', energies)
        np.savetxt(savedir + '/iterations.txt', iterations)

    # Examples A/1 and B/2 in Gradient Damage Models
    U_t = L*t

    # Example A/1 in Gradient Damage Models
    if case == "A":
        if U_t == 0:
            alpha_t = -math.inf
        else:
            alpha_t = 1-(U_e/U_t)**2

        alpha_t = max(0, alpha_t)

        if U_t <= U_e:
            sigma_t = sigma_e*(U_t/U_e)
        else:
            sigma_t = sigma_e*(U_e/U_t)**3

        # Save dimensionalized data
        dim_data[i_t] = np.array([U_t/U_e, alpha_t, sigma_t/sigma_e])

    # Example B/1 in Gradient Damage Models
    else:
        alpha_t = U_t**2/(3*U_m**2 + U_t**2)
        sigma_t = (9*E_0/L)*(U_t*U_m**4)/((3*U_m**2 + U_t**2)**2)

        # Save dimensionalized data
        dim_data[i_t] = np.array([U_t/U_m, alpha_t, sigma_t/sigma_m])

    # Save non-dimensionalized data
    nondim_data[i_t] = np.array([U_t, alpha_t, sigma_t])

# Plotting
# ----------------------------------------------------------------------------
if MPI.rank(MPI.comm_world) == 0:
    plt.figure(1)
    p1, = plt.plot(energies[:, 0], energies[:, 1])
    p2, = plt.plot(energies[:, 0], energies[:, 2])
    p3, = plt.plot(energies[:, 0], energies[:, 3])
    plt.legend([p1, p2, p3], ["Elastic", "Dissipated", "Total"], loc="best", frameon=False)
    plt.xlabel('Displacement')
    plt.ylabel('Energies')
    plt.savefig(savedir + '/energies.pdf', transparent=True)
    plt.close()

    plt.figure(2)
    plt.plot(nondim_data[:, 0], nondim_data[:, 1], 'k*-')
    plt.xlabel('Strain (U_t)')
    plt.ylabel('Damage (alpha_t)')
    plt.title('Damage vs Strain ' + filename)
    plt.savefig(savedir + '/damage_strain_nondim.pdf', transparent=True)
    plt.close()

    plt.figure(3)
    plt.plot(dim_data[:, 0], dim_data[:, 1], 'k*-')
    plt.xlabel('Strain (U_t/U_e)')
    plt.ylabel('Damage (alpha_t)')
    plt.title('Damage vs Strain ' + filename)
    plt.savefig(savedir + '/damage_strain_dim.pdf', transparent=True)
    plt.close()

    plt.figure(4)
    plt.plot(nondim_data[:, 0], nondim_data[:, 2], 'k*-')
    plt.xlabel('Strain (U_t/U_m)')
    plt.ylabel('Stress (sigma_t/sigma_m)')
    plt.title('Stress vs Strain ' + filename)
    plt.savefig(savedir + '/stress_strain_nondim.pdf', transparent=True)
    plt.close()

    plt.figure(3)
    plt.plot(dim_data[:, 0], dim_data[:, 2], 'k*-')
    plt.xlabel('Strain (U_t)')
    plt.ylabel('Stress (sigma_t)')
    plt.title('Stress vs Strain ' + filename)
    plt.savefig(savedir + '/stress_strain_dim.pdf', transparent=True)
    plt.close()
