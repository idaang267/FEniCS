#
# =============================================================================
# FEniCS code for Variational Fracture Mechanics
# =============================================================================
#
# A static solution of the variational fracture mechanics problems
# using the regularization two-fold anisotropic damage model
#
# author: Bin Li (bl736@cornell.edu), Corrado Maurini (corrado.maurini@upmc.fr)
#
# date: 10/10/2017
# ----------------


# ----------------------------------------------------------------------------
from __future__ import division
from dolfin import *
from mshr import *
from ufl import replace

import argparse
import math
import os
import shutil
import sympy
import sys
import numpy as np
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Parameters for DOLFIN and SOLVER
# ----------------------------------------------------------------------------
#set_log_level(LogLevel.INFO)  # log level
set_log_level(LogLevel.WARNING)
# set some dolfin specific parameters
info(parameters,True)
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"

# -----------------------------------------------------------------------------
# parameters of the solvers
"""
solver_u_parameters = {"nonlinear_solver": "newton",
                       "newton_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 100,
                                          "absolute_tolerance": 1e-8,
                                          "relative_tolerance": 1e-6,
                                          "report": True,
                                          "error_on_nonconvergence": True}}
"""
solver_u_parameters = {"linear_solver": "mumps",
                       "symmetric": True,
                       "preconditioner": "hypre_amg",
                       "krylov_solver": {
                           "report": False,
                           "monitor_convergence": False,
                           "relative_tolerance": 1e-10,
                           "absolute_tolerance": 1e-10}}
# parameters of the PETSc/Tao solver used for the alpha-problem
tao_solver_parameters = {"maximum_iterations": 200,
                         "report": False,
                         "line_search": "more-thuente",
                         "linear_solver": "mumps",
                         "method": "tron",
                         "gradient_absolute_tol": 1e-8,
                         "gradient_relative_tol": 1e-8,
                         "error_on_nonconvergence": False}

# -----------------------------------------------------------------------------
# set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("meshsize", 200)
userpar.add("load_min", 0.)
userpar.add("load_max", 4.0)
userpar.add("load_steps", 201)
userpar.parse()

# ----------------------------------------------------------------------------
# Parameters for ANISOTROPIC surface energy and materials
# ----------------------------------------------------------------------------

# Material constant
E       = Constant(1.0)
nu      = Constant(0.3)
Gc      = Constant(1.0)
k_ell   = Constant(1.e-6)  # residual stiffness

# Loading Parameters
ut      = 1.0   # reference value for the loading (imposed displacement)
w1      = 1.0
E_0     = 1.0

# Numerical parameters of the alternate minimization
maxiteration = 2000
AM_tolerance = 1e-4

# Geometry paramaters
L = 1.0
H = 0.1
N         = userpar["meshsize"]
hsize     = float(L/N)
ell       = Constant(6.0*hsize) # damage paramaters

# must edit between ExampleA and ExampleB
case = "B"
modelname = "traction_bar"
meshname  = modelname+"_mesh.xdmf"
simulation_params = "h_%.4f" % (hsize)
savedir   = "output" + case + "/" + simulation_params + "/"
filename = "Example" + case

if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

mesh     = RectangleMesh(Point(0., -0.5*H), Point(L, 0.5*H), int(N), int(float(H/hsize)), "right/left")
#geometry = Rectangle(Point(0., -0.5*H), Point(L, 0.5*H))
#mesh     = generate_mesh(geometry, int(N), 'cgal')
geo_mesh  = XDMFFile(MPI.comm_world, savedir+meshname)
geo_mesh.write(mesh)

mesh.init()
ndim = mesh.geometry().dim()  # get number of space dimensions
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# ----------------------------------------------------------------------------
# Strain and stress and Constitutive functions of the damage model
# ----------------------------------------------------------------------------
# Strain and stress
def eps(v):
    return sym(grad(v))

def sigma_0(v):
    mu = E/(2.0*(1.0+nu))
    lmbda = E*nu/(1.0-nu**2)  # plane stress
    return 2.0*mu*(eps(v))+lmbda*tr(eps(v))*Identity(ndim)

# Constitutive functions of the damage model
def w(alpha):
    w1 = 1.0
    if case == "A":
        return w1*alpha
    else:
        return w1*alpha**2

def a(alpha):
    return (1.0-alpha)**2

# ----------------------------------------------------------------------------
# Define boundary sets for boundary conditions
# Impose the displacements field
# ----------------------------------------------------------------------------
def left_boundary(x, on_boundary):
    return on_boundary and near(x[0], 0.0, 0.1 * hsize)

def right_boundary(x, on_boundary):
    return on_boundary and near(x[0], L, 0.1 * hsize)

## when using "pointwise", the boolean argument on_boundary
## in SubDomain::inside will always be false
def left_pinpoints(x, on_boundary):
    return  near(x[0], 0.0, 0.1 * hsize) and near(x[1], 0.0, 0.1 * hsize)
def right_pinpoints(x, on_boundary):
    return  near(x[0], L, 0.1 * hsize) and near(x[1], 0.0, 0.1 * hsize)

# ----------------------------------------------------------------------------
# Variational formulation
# ----------------------------------------------------------------------------
# Create function space for 2D elasticity + Damage
V_u       = VectorFunctionSpace(mesh, "Lagrange", 1)
V_alpha   = FunctionSpace(mesh, "Lagrange", 1)
V_sigma   = TensorFunctionSpace(mesh, "Lagrange", 1)
V_epsilon = TensorFunctionSpace(mesh, "Lagrange", 1)

# Define the function, test and trial fields
u      = Function(V_u, name="Displacement")
du     = TrialFunction(V_u)
v      = TestFunction(V_u)
alpha  = Function(V_alpha, name="Damage")
dalpha = TrialFunction(V_alpha)
beta   = TestFunction(V_alpha)

# --------------------------------------------------------------------
# Dirichlet boundary condition
# Impose the displacements field given by asymptotic expansion of crack tip
# --------------------------------------------------------------------
u_UL = Expression("0.0", degree=0)
#u_UL = Expression(["0.0", "0.0"], degree=0)
u_UR = Expression("t", t=0.0, degree=0) # slide Dirichlet BCs

# bc - u (imposed displacement) # slide Dirichlet boundary condition
Gamma_u_0 = DirichletBC(V_u.sub(0), u_UL, left_boundary)
Gamma_u_1 = DirichletBC(V_u.sub(1), u_UL, left_pinpoints, method='pointwise')
Gamma_u_2 = DirichletBC(V_u.sub(0), u_UR, right_boundary)
#Gamma_u_3 = DirichletBC(V_u.sub(1), u_UL, right_pinpoints, method='pointwise')
bc_u = [Gamma_u_0, Gamma_u_1, Gamma_u_2]

# bc - alpha (zero damage)
Gamma_alpha_0 = DirichletBC(V_alpha, 0.0, left_boundary)
Gamma_alpha_1 = DirichletBC(V_alpha, 0.0, right_boundary)
bc_alpha      = [Gamma_alpha_0,Gamma_alpha_1]

# --------------------------------------------------------------------
# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Fenics forms for the energies
def sigma(u, alpha):
    return (a(alpha)+k_ell)*sigma_0(u)

body_force        = Constant((0., 0.))
elastic_energy    = 0.5*inner(sigma(u, alpha), eps(u))*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy-external_work

# Weak form of elasticity problem
"""
# Writing tangent problems in term of test and trial functions for matrix assembly
E_du = derivative(E_u, u, du)

# Variational problem for the displacement
problem_u = NonlinearVariationalProblem(E_u, u, bc_u, J=E_du)
# Set up the solvers
solver_u  = NonlinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
# info(solver_u.parameters, True)
"""
E_u = derivative(elastic_potential, u, v)
# Writing tangent problems in term of test and trial functions for matrix assembly
E_du = replace(E_u, {u: du})
# Variational problem for the displacement
problem_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), u, bc_u)
# Set up the solvers
solver_u  = LinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
# info(solver_u.parameters, True)

# --------------------------------------------------------------------
# Define the energy functional of damage problem
# --------------------------------------------------------------------
alpha_0 = interpolate(Expression("0.", degree=0), V_alpha)  # initial (known) alpha
z = sympy.Symbol("z", positive=True)
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*inner(grad(alpha), grad(alpha)))*dx
damage_functional = elastic_potential + dissipated_energy

# Compute directional derivative about alpha in the direction of beta (Gradient)
E_alpha       = derivative(damage_functional, alpha, beta)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
E_alpha_alpha = derivative(E_alpha, alpha, dalpha)
# --------------------------------------------------------------------
# Implement the box constraints for damage field
# --------------------------------------------------------------------
# Variational problem for the damage (non-linear to use variational inequality solvers of petsc)
# Define the minimisation problem by using OptimisationProblem class
class DamageProblem(OptimisationProblem):

    def __init__(self,f,gradf,alpha,J,bcs):
        OptimisationProblem.__init__(self)
        self.total_energy = f
        self.Dalpha_total_energy = gradf
        self.J_alpha = J
        self.alpha = alpha
        self.bc_alpha = bcs

    def f(self, x):
        self.alpha.vector()[:] = x
        return assemble(self.total_energy)

    def F(self, b, x):
        self.alpha.vector()[:] = x
        assemble(self.Dalpha_total_energy, b)
        for bc in self.bc_alpha:
            bc.apply(b)

    def J(self, A, x):
        self.alpha.vector()[:] = x
        assemble(self.J_alpha, A)
        for bc in self.bc_alpha:
            bc.apply(A)

damage_problem = DamageProblem(damage_functional,E_alpha,alpha,E_alpha_alpha,bc_alpha)

# Set up the solvers
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
alpha_lb = interpolate(Expression("0.", degree=0), V_alpha)  # lower bound, set to 0
alpha_ub = interpolate(Expression("1.", degree=0), V_alpha)  # upper bound, set to 1

# loading and initialization of vectors to store time datas
load_multipliers      = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies              = np.zeros((len(load_multipliers), 4))
iterations            = np.zeros((len(load_multipliers), 2))

save_alpha_t           = np.zeros((len(load_multipliers), 2))
save_alpha_t_nondim    = np.zeros((len(load_multipliers), 2))
save_sigma_t           = np.zeros((len(load_multipliers), 2))
save_sigma_t_nondim    = np.zeros((len(load_multipliers), 2))

# set the saved data file name
results      = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
results.parameters["rewrite_function_mesh"]          = False
results.parameters["flush_output"]                   = True
results.parameters["functions_share_mesh"]           = True
# write the parameters to file
File(savedir+"/parameters.xml") << userpar

# Constants for Examples A/1 and B/2 for Gradient Damage Models
if case == "A":
    sigma_e = sqrt(w1*E_0)
    U_e = (sigma_e/E_0)*L
else:
    sigma_m = (3*sqrt(3))/(8*sqrt(2))*sqrt(w1*E_0)
    U_m = (16/9)*(sigma_m/E_0)*L

# ----------------------------------------------------------------------------
# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    u_UR.t = t * ut # here, ut = L
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))
    # Alternate Mininimization
    # Initialization
    iteration = 1
    err_alpha = 1.0
    # Iterations
    while err_alpha > AM_tolerance and iteration < maxiteration:
        # solve elastic problem
        solver_u.solve()
        # solve damage problem with box constraint
        solver_alpha.solve(damage_problem, alpha.vector(), alpha_lb.vector(), alpha_ub.vector())
        # test error
        alpha_error = alpha.vector() - alpha_0.vector()
        err_alpha = alpha_error.norm('linf')
        # monitor the results
        if MPI.rank(MPI.comm_world) == 0:
          print ("AM Iteration: {0:3d},  alpha_error: {1:>14.8f}".format(iteration, err_alpha))

        # update iteration
        alpha_0.assign(alpha)
        iteration = iteration + 1
    # updating the lower bound to account for the irreversibility
    alpha_lb.vector()[:] = alpha.vector()
    alpha.rename("Damage", "alpha")
    u.rename("Displacement", "u")

    sigma_val = project(sigma(u, alpha), V_sigma)
    sigma_val.rename("Stress", "sigma_val")

    epsilon_val = project(eps(u), V_epsilon)
    epsilon_val.rename("Strain", "epsilon_val")

    # Dump solution to file
    results.write(alpha, t)
    results.write(u, t)
    results.write(sigma_val, t)
    results.write(epsilon_val, t)

    # ----------------------------------------
    # Some post-processing
    # ----------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value = assemble(elastic_energy)
    surface_energy_value = assemble(dissipated_energy)
    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, elastic_energy_value+surface_energy_value])

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

        save_alpha_t[i_t]        = np.array([U_t, alpha_t])
        save_alpha_t_nondim[i_t]  = np.array([U_t/U_e, alpha_t])
        save_sigma_t[i_t]        = np.array([U_t, sigma_t])
        save_sigma_t_nondim[i_t] = np.array([U_t/U_e, sigma_t/sigma_e])

    # Example B/1 in Gradient Damage Models
    else:
        alpha_t = U_t**2/(3*U_m**2 + U_t**2)
        sigma_t = (9*E/L)*(U_t*U_m**4)/((U_m**2 + U_t**2)**2)

        save_alpha_t[i_t]        = np.array([U_t, alpha_t])
        save_alpha_t_nondim[i_t] = np.array([U_t/U_m, alpha_t])
        save_sigma_t[i_t]        = np.array([U_t, sigma_t])
        save_sigma_t_nondim[i_t] = np.array([U_t/U_m, sigma_t/sigma_m])

# ----------------------------------------------------------------------------
# Plot energy and stresses
if MPI.rank(MPI.comm_world) == 0:
    plt.figure(1)
    p1, = plt.plot(energies[:, 0], energies[:, 1])
    p2, = plt.plot(energies[:, 0], energies[:, 2])
    p3, = plt.plot(energies[:, 0], energies[:, 3])
    plt.legend((p1, p2, p3), ('Elastic', 'Dissipated', 'Total'), loc="best", frameon=False)
    plt.xlabel('Displacement')
    plt.ylabel('Energies')
    plt.title('Energy vs Displacement ' + filename)
    plt.savefig(savedir + '/energies.pdf', transparent=True)
    plt.close()

    plt.figure(2)
    plt.plot(save_alpha_t_nondim[:, 0], save_alpha_t_nondim[:, 1], 'k*-')
    plt.xlabel('Strain (U_t/U_e)')
    plt.ylabel('Damage (alpha_t)')
    plt.title('Damage vs Strain ' + filename)
    plt.savefig(savedir + '/damage_strain_nondim.pdf', transparent=True)
    plt.close()

    plt.figure(3)
    plt.plot(save_sigma_t[:, 0], save_sigma_t[:, 1], 'k*-')
    plt.xlabel('Strain (U_t)')
    plt.ylabel('Stress (sigma_t)')
    plt.title('Stress vs Strain ' + filename)
    plt.savefig(savedir + '/stress_strain.pdf', transparent=True)
    plt.close()

    plt.figure(4)
    plt.plot(save_sigma_t_nondim[:, 0], save_sigma_t_nondim[:, 1], 'k*-')
    plt.xlabel('Strain (U_t/U_m)')
    plt.ylabel('Stress (sigma_t/sigma_m)')
    plt.title('Stress vs Strain ' + filename)
    plt.savefig(savedir + '/stress_strain_nondim.pdf', transparent=True)
    plt.close()
