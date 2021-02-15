# -------------------------------------------
# FEniCS code  Variational Fracture Mechanics
################################################################################
#                                                                              #
# A stabilized mixed finite element method for gradient damage models of 	   #
# fracture in incompressible hyperelastic materials                            #
# Author: Bin Li                                                               #
# Email: bl736@cornell.edu                                                     #
# date: 10/01/2018                                                             #
#                                                                              #
################################################################################
# e.g. python3 traction-stabilized.py --meshsize 100						   #
################################################################################

# ----------------------------------------------------------------------------
from __future__ import division
from dolfin import *
from mshr import *
from scipy import optimize

import argparse
import math
import os
import shutil
import sympy
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Parameters for DOLFIN and SOLVER
# ----------------------------------------------------------------------------
set_log_level(LogLevel.WARNING)  # 20 Information of general interest

# Set some dolfin specific parameters
# ----------------------------------------------------------------------------
parameters["form_compiler"]["representation"]="uflacs"
parameters["form_compiler"]["optimize"]=True
parameters["form_compiler"]["cpp_optimize"]=True
parameters["form_compiler"]["quadrature_degree"]=2
info(parameters,True)

# Parameters of the solvers for displacement and damage (alpha-problem)
# -----------------------------------------------------------------------------
# Parameters of the nonlinear SNES solver used for the displacement u-problem
solver_up_parameters  = {"nonlinear_solver": "snes",
                         "symmetric": True,
                         "snes_solver": {"linear_solver": "mumps",
                                         "method" : "newtontr",
                                         "line_search": "cp",
                                         "preconditioner" : "hypre_amg",
                                         "maximum_iterations": 200,
                                         "absolute_tolerance": 1e-10,
                                         "relative_tolerance": 1e-10,
                                         "solution_tolerance": 1e-10,
                                         "report": True,
                                         "error_on_nonconvergence": False}}

# Parameters of the PETSc/Tao solver used for the alpha-problem
tao_solver_parameters = {"maximum_iterations": 100,
                         "report": False,
                         "line_search": "more-thuente",
                         "linear_solver": "cg",
                         "preconditioner" : "hypre_amg",
                         "method": "tron",
                         "gradient_absolute_tol": 1e-8,
                         "gradient_relative_tol": 1e-8,
                         "error_on_nonconvergence": True}

# --------------------------------------------------------------------
# Implement the box constraints for damage field
# --------------------------------------------------------------------
# Variational problem for the damage problem
# (non-linear - use variational inequality solvers in PETSc)

# Define the minimisation problem by using OptimisationProblem class
class DamageProblem(OptimisationProblem):
    def __init__(self):
        OptimisationProblem.__init__(self)
        self.total_energy = damage_functional
        self.Dalpha_total_energy = E_alpha
        self.J_alpha = E_alpha_alpha
        self.alpha = alpha
        self.bc_alpha = bc_alpha
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

# Set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("mu", 1)           # Shear modulus
userpar.add("kappa",10e2)      # Bulk modulus
userpar.add("Gc", 2.4E6)       # Fracture toughness (2.4E3)
userpar.add("k_ell", 5.e-5)    # Residual stiffness
userpar.add("meshsizeX", 100)
userpar.add("meshsizeY", 50)
userpar.add("meshsizeZ", 50)
userpar.add("load_min", 0.)
userpar.add("load_max", 0.01)
userpar.add("load_steps", 10)
userpar.add("ell_multi", 5)
# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
# Geometry parameters
L, H, W = 5, 1, 1            # Length (x), height (y), and width (x-direction)
Nx   = userpar["meshsizeX"]
Ny   = userpar["meshsizeY"]
Nz   = userpar["meshsizeZ"]
hsize = float(H/Ny)            # Geometry based definition for regularization
S = userpar["load_steps"]

# Shear modulus and Poisson's ratio
mu    = userpar["mu"]           # Shear modulus
kappa = userpar["kappa"]        # Bulk Modulus

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]
ell_multi = userpar["ell_multi"]
# Damage regularization parameter - internal length scale used for tuning Gc
ell   = Constant(ell_multi*hsize)

# Naming parameters for saving output
modelname = "3D-stabilized"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "L_%.0f_Nx_%.0f_H_%.0f_Ny_%.0f_W_%.0f_Nz_%.0f_lx_%.0f_d_%.2f_k_%.0f" % (L, Nx, H, Ny, W, Nz, ell_multi, userpar["load_max"], kappa)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

# For parallel processing - write one directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(L, H, W), Nx, Ny, 10)
# mesh = Mesh("ShearTestRef.xml")
geo_mesh = XDMFFile(MPI.comm_world, savedir + meshname)
geo_mesh.write(mesh)

# Obtain number of space dimensions
mesh.init()
ndim = mesh.geometry().dim()
# Structure used for one printout of the statement
if MPI.rank(MPI.comm_world) == 0:
    print ("Mesh Dimension: {0:2d}".format(ndim))

# Characteristic element length - used for stabilization
h = CellDiameter(mesh)

# Reference value for the loading (imposed displacement)
ut = 1.0

# Numerical parameters of the alternate minimization scheme
maxiteration = 2000         # Sets a limit on number of iterations
AM_tolerance = 1e-4

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class left_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, 0.1 * hsize)

class right_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, 0.1 * hsize)

class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0, 0.1 * hsize)

class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, 0.1 * hsize) #and between(x[0], (5, L))

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L/2, 0.1*hsize) and near(x[1], H/2, 0.1*hsize) and near(x[2], W/2, 0.1*hsize)

# Convert all boundary classes for visualization
left_boundary = left_boundary()
right_boundary = right_boundary()
top_boundary = top_boundary()
bot_boundary = bot_boundary()
pin_point = pin_point()

lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 3)

# show lines of interest
lines.set_all(0)
bot_boundary.mark(lines, 1)
top_boundary.mark(lines, 1)
file_results = XDMFFile(savedir + "/" + "lines.xdmf")
file_results.write(lines)

# Show points of interest
points.set_all(0)
pin_point.mark(points, 1)
file_results = XDMFFile(savedir + "/" + "points.xdmf")
file_results.write(points)

# Variational formulation
# ----------------------------------------------------------------------------
# Create function space for elasticity + damage
P1  = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
P2  = VectorElement("Lagrange", mesh.ufl_cell(), 1)
# Stabilized mixed FEM for incompressible elasticity
TH  = MixedElement([P2,P1])
# Define function spaces for displacement and pressure in V_u and damage in V_alpha
V_u = FunctionSpace(mesh, TH)
V_alpha = FunctionSpace(mesh, "Lagrange", 1)

# Define the function, test and trial fields for elasticity problem
w_p    = Function(V_u)
u_p    = TrialFunction(V_u)
v_q    = TestFunction(V_u)
(u, p) = split(w_p)     # Displacement and pressure (u and p)
(v, q) = split(v_q)     # Test functions for u and p
# Define the function, test and trial fields for damage problem
alpha  = Function(V_alpha, name = "Damage")
dalpha = TrialFunction(V_alpha)
beta   = TestFunction(V_alpha)

# Dirichlet boundary condition
# --------------------------------------------------------------------
u00 = Constant((0.0))
u0 = Expression(["0.0", "0.0", "0.0"], degree=0)
u1 = Expression("t", t= 0.0, degree=0)
u2 = Expression("-t", t= 0.0, degree=0)
# bc - u (imposed displacement)
# roller boundary
bc_u0 = DirichletBC(V_u.sub(0).sub(0), u00, right_boundary)
# top and bottom boundaries are subjected to a displacement in the y direction
bc_u1 = DirichletBC(V_u.sub(0).sub(1), u1, top_boundary)
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u2, bot_boundary)
# Pinned line in the center
bc_u5 = DirichletBC(V_u.sub(0), u0, pin_point, method="pointwise")
# Combine boundary conditions
bc_u = [bc_u1, bc_u2, bc_u5]

# bc - alpha (zero damage)
# No damage to the boundaries - damage does not initiate from constrained edges
bc_alpha0 = DirichletBC(V_alpha, 0.0, bot_boundary)
bc_alpha1 = DirichletBC(V_alpha, 0.0, top_boundary)
bc_alpha = [bc_alpha0, bc_alpha1]

# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Define some parameters for the eigenvalues
d_par = tr(C)/3.
e_par = sqrt(tr(C-d_par*I)**2/6.)
f_par = (1./e_par)*(C-d_par*I)
g_par = det(f_par)/2

# Define the eigenvalues of C (principal stretches)
lambda_1 = sqrt(d_par + 2.*e_par*cos(acos(g_par)/3.))
lambda_3 = sqrt(d_par + 2.*e_par*cos((acos(g_par)/3.) + 2.*pi/3.))
lambda_2 = sqrt(3.*d_par - lambda_1**2 - lambda_3**3)

# m_par = tr(C)/3
# q_par = det(C-m_par*I)/2
# p_par = tr(C-m_par*I)**2/6
#
# phi = 1/3*atan(sqrt(p_par**3 - q_par**2)/q_par)
# lambda_1 = m_par + 2*sqrt(p_par)*cos(phi)
# lambda_2 = m_par - sqrt(p_par)*(cos(phi) + sqrt(3)*sin(phi))
# lambda_3 = m_par - sqrt(p_par)*(cos(phi) - sqrt(3)*sin(phi))

# Define the energy functional of the elasticity problem
# --------------------------------------------------------------------

# Functions of the damage model
# ----------------------------------------------------------------------------
def hs(x):# Heaviside function
    return (x + abs(x))/(2.*x)

hs_p_l1 = conditional(gt(lambda_1,1), lambda_1, 1)
hs_p_l2 = conditional(gt(lambda_2,1), lambda_2, 1)
hs_p_l3 = conditional(gt(lambda_3,1), lambda_3, 1)

hs_p_J = conditional(gt(J,1), J, 1)

hs_n_l1 = conditional(lt(lambda_1,1), lambda_1, 1)
hs_n_l2 = conditional(lt(lambda_2,1), lambda_2, 1)
hs_n_l3 = conditional(lt(lambda_3,1), lambda_3, 1)

hs_n_J = conditional(lt(J,1), J, 1)

def w(alpha):
    return alpha

def a(alpha):
    return (1.0-alpha)**2

def b(alpha):
    return (1.0-alpha)**3

# Zero body force
body_force = Constant((0., 0., 0.))
# Non-dimension non-negative stability parameter
varpi_ = 1.0
# Eq 19 in Klaas
varpi  = project(varpi_*h**2/(2.0*mu), FunctionSpace(mesh,'DG',0))
# Elastic energy
W_act = (a(alpha)+k_ell)*(mu/2.)*(hs_p_l1**2-1.-2.*ln(hs_p_l1))*dx \
      + (a(alpha)+k_ell)*(mu/2.)*(hs_p_l2**2-1.-2.*ln(hs_p_l2))*dx \
      + (a(alpha)+k_ell)*(mu/2.)*(hs_p_l3**2-1.-2.*ln(hs_p_l3))*dx \
      + a(alpha)**3*(kappa/2.)*(hs_p_J-1.)**2*dx
W_pas = (a(alpha)+k_ell)*(mu/2.)*(hs_n_l1**2-1.-2.*ln(hs_n_l1))*dx \
      + (a(alpha)+k_ell)*(mu/2.)*(hs_n_l2**2-1.-2.*ln(hs_n_l2))*dx \
      + (a(alpha)+k_ell)*(mu/2.)*(hs_n_l3**2-1.-2.*ln(hs_n_l3))*dx \
      + a(alpha)**3*(kappa/2.)*(hs_n_J-1.)**2*dx
# additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy    = W_act + W_pas
# elastic_energy    = (a(alpha)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
#                     - b(alpha)*p*(J-1.0)*dx - 1/(2*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Define the stabilization term
# Compute directional derivative about w_p in the direction of v (Gradient)
F_u = derivative(elastic_potential, w_p, v_q) #\ - varpi*b(alpha)*J*inner(inv(C), outer(grad(p),grad(q)))*dx
# Compute directional derivative about w_p in the direction of u_p (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem to solve for displacement and pressure
problem_up = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solver for displacement and pressure
solver_up  = NonlinearVariationalSolver(problem_up)
solver_up.parameters.update(solver_up_parameters)
# info(solver_up.parameters, True) # uncomment to see available parameters

# Define the energy functional of the damage problem
# --------------------------------------------------------------------
# Initializing known alpha
alpha_0 = interpolate(Expression("0.", degree=0), V_alpha)
# Define the specific energy dissipation per unit volume
z = sympy.Symbol("z", positive=True)
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
# Define the phase-field fracture term of the damage functional
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*dot(grad(alpha), grad(alpha)))*dx
damage_functional = elastic_potential + dissipated_energy

# Compute directional derivative about alpha in the direction of beta (Gradient)
E_alpha = derivative(damage_functional, alpha, beta)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
E_alpha_alpha = derivative(E_alpha, alpha, dalpha)

# Set the lower and upper bound of the damage variable (0-1)
# alpha_lb = interpolate(Expression("0.", degree=0), V_alpha)
alpha_lb = interpolate(Expression("x[0]>=0 & x[0]<=2.5 & near(x[1], 0.5, 0.1 * hsize) ? 1.0 : 0.0", \
                       hsize = hsize, degree=0), V_alpha)
alpha_ub = interpolate(Expression("1.", degree=0), V_alpha)

# Set up the solvers
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
# info(solver_alpha.parameters, True) # uncomment to see available parameters

# Loading and initialization of vectors to store data of interest
load_multipliers = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies         = np.zeros((len(load_multipliers), 5))
iterations       = np.zeros((len(load_multipliers), 2))

# Split into displacement and pressure
(u, p)   = w_p.split()
# Data file name
file_tot = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
# Saves the file in case of interruption
file_tot.parameters["rewrite_function_mesh"] = False
file_tot.parameters["functions_share_mesh"]  = True
file_tot.parameters["flush_output"]          = True
# Write the parameters to file
File(savedir+"/parameters.xml") << userpar

timer0 = time.process_time()

# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    # Update the displacement with each iteration
    u1.t = t * ut
    u2.t = t * ut
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))

    # Alternate Mininimization scheme
    # -------------------------------------------------------------------------
    # Solve for u holding alpha constant then solve for alpha holding u constant
    iteration = 1           # Initialization of iteration loop
    err_alpha = 1.0         # Initialization for condition for iteration
    # Conditions for iteration
    while err_alpha > AM_tolerance and iteration < maxiteration:
        # Solve elastic problem
        solver_up.solve()
        # Solve damage problem with box constraint
        solver_alpha.solve(DamageProblem(), alpha.vector(), alpha_lb.vector(), alpha_ub.vector())
        # Update the alpha condition for iteration by calculating the alpha error norm
        alpha_error = alpha.vector() - alpha_0.vector()
        err_alpha = alpha_error.norm('linf')    # Row-wise norm
        # Printouts to monitor the results and number of iterations
        volume_ratio = assemble(J/(L*H)*dx)
        if MPI.rank(MPI.comm_world) == 0:
            print ("AM Iteration: {0:3d},  alpha_error: {1:>14.8f}".format(iteration, err_alpha))
            print("\nVolume Ratio: [{}]".format(volume_ratio))
        # Update variables for next iteration
        alpha_0.assign(alpha)
        iteration = iteration + 1

    # Updating the lower bound to account for the irreversibility of damage
    alpha_lb.vector()[:] = alpha.vector()

    # Rename for paraview
    alpha.rename("Damage", "alpha")
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")

    # Write solution to file
    file_tot.write(alpha, t)
    file_tot.write(u, t)
    file_tot.write(p, t)

    # Post-processing
    # --------------------------------------------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value = assemble(elastic_energy)
    surface_energy_value = assemble(dissipated_energy)

    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, \
                              elastic_energy_value+surface_energy_value, volume_ratio])

    if MPI.rank(MPI.comm_world) == 0:
        print("\nEnd of timestep {0:3d} with load multiplier {1:4f}".format(i_t, t))
        print("\nElastic and Surface Energies: [{0:6f},{1:6f}]".format(elastic_energy_value, surface_energy_value))
        print("\nElastic and Surface Energies: [{},{}]".format(elastic_energy_value, surface_energy_value))
        print("\nVolume Ratio: [{}]".format(volume_ratio))
        print("-----------------------------------------")
        # Save some global quantities as a function of the time
        np.savetxt(savedir + '/stabilized-energies.txt', energies)
        np.savetxt(savedir + '/stabilized-iterations.txt', iterations)
# ----------------------------------------------------------------------------

print("elapsed CPU time: ", (time.process_time() - timer0))


# # Plot energy and stresses
# if MPI.rank(MPI.comm_world) == 0:
#     p1, = plt.plot(energies[slice(None), 0], energies[slice(None), 1])
#     p2, = plt.plot(energies[slice(None), 0], energies[slice(None), 2])
#     p3, = plt.plot(energies[slice(None), 0], energies[slice(None), 3])
#     plt.legend([p1, p2, p3], ["Elastic", "Dissipated", "Total"], loc="best", frameon=False)
#     plt.xlabel('Displacement')
#     plt.ylabel('Energies')
#     plt.title('stabilized FEM')
#     plt.savefig(savedir + '/stabilized-energies.pdf', transparent=True)
#     plt.close()
#     p4, = plt.plot(energies[slice(None), 0], energies[slice(None), 4])
#     plt.xlabel('Displacement')
#     plt.ylabel('Volume ratio')
#     plt.title('stabilized FEM')
#     plt.savefig(savedir + '/stabilized-volume-ratio.pdf', transparent=True)
#     plt.close()
