# FEniCS code  Variational Fracture Mechanics
#################################################################################################################
#                                                                                                               #
# A Taylor-Hood finite element method for gradient damage models of                                             #
# fracture in incompressible hyperelastic materials
# Addition: Hybrid formulation
# author: Bin Li                                                                                                #
# Email: bl736@cornell.edu                                                                                      #
# date: 2/21/2021                                                                                              #
#                                                                                                               #
#################################################################################################################
# e.g. python3 traction-neo-Hookean.py --meshsize 100														                              	#
#################################################################################################################

# Import Modules
# ----------------------------------------------------------------------------
from __future__ import division
from dolfin import *
from mshr import *
from scipy import optimize
from ufl import eq

import argparse
import math
import os
import shutil
import sympy
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

# Parameters for DOLFIN and SOLVER
# ----------------------------------------------------------------------------
set_log_level(LogLevel.WARNING)  # 20, // information of general interest

# set some dolfin specific parameters
# ----------------------------------------------------------------------------
parameters["form_compiler"]["representation"]="uflacs"
parameters["form_compiler"]["optimize"]=True
parameters["form_compiler"]["cpp_optimize"]=True
parameters["form_compiler"]["quadrature_degree"]=2
info(parameters,True)

# Parameters of the solvers for displacement and damage (alpha-problem)
# -----------------------------------------------------------------------------
# Parameters of the nonlinear SNES solver used for the displacement u-problem
solver_u_parameters   = {"nonlinear_solver": "snes",
                         "symmetric": True,
                         "snes_solver": {"linear_solver": "mumps", #"mumps",
                                         "method" : "newtontr", #"newtontr", "newtonls"
                                         "line_search": "cp", # "basic", "bt", "l2", "cp", "nleqerr"
                                         "preconditioner" : "hypre_amg",
                                         "maximum_iterations": 200,
                                         "absolute_tolerance": 1e-10,
                                         "relative_tolerance": 1e-10,
                                         "solution_tolerance": 1e-10,
                                         "report": True,
                                         "error_on_nonconvergence": True}}

# parameters of the PETSc/Tao solver used for the alpha-problem
tao_solver_parameters = {"maximum_iterations": 100,
                         "line_search": "more-thuente",
                         "linear_solver": "cg",
                         "preconditioner" : "hypre_amg",
                         "method": "tron",
                         "gradient_absolute_tol": 1e-8,
                         "gradient_relative_tol": 1e-8,
                         "report": False,
                         "error_on_nonconvergence": True}

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

# set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("mu", 0.2)         # Shear Modulus
userpar.add("kappa", 200)      # Bulk Modulus
userpar.add("Gc", 1)           # fracture toughness
userpar.add("k_ell", 5.e-5)    # residual stiffness
userpar.add("load_min", 0.00)
userpar.add("load_max", 3.5)
userpar.add("load_steps", 400)
userpar.add("ell_multi", 5)
# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
# Geometry parameters
W, H, T = 1.0, 1.5, 0.02
hsize   = 0.01         # Element size float(L/N)

# Material parameters
mu    = userpar["mu"]
kappa = userpar["kappa"]

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]
ell_multi = userpar["ell_multi"]
# Damage regularization parameter - internal length scale used for tuning Gc
ell = Constant(ell_multi*hsize)

# Number of steps
S = userpar["load_steps"]
# Numerical parameters of the alternate minimization scheme
maxiteration = 4000
AM_tolerance = 1e-4

# Naming parameters for saving output
modelname = "3D-Stabilized"
meshname  = modelname+"-mesh.xdmf"
simulation_params = "Form3_W_%.1f_H_%.1f_T_%.3f_h_%.3f_d_%.2f_S_%.0f" % (W, H, T, hsize, userpar["load_max"],S)
savedir   = "output/"+modelname+"/"+simulation_params+"/"

# For parallel processing - write one directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Different mesh generation methods
# ----------------------------------------------------------------------------
# Mesh generation for geo file
mesh = Mesh("3DEdgeCrack.xml")
geo_mesh  = XDMFFile(MPI.comm_world, savedir+meshname)
geo_mesh.write(mesh)

# # Mesh generation for structured mesh
# mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(L, H, W), N, int(N*H/L), int(N*W/L))

# # Mesh generation from xdmf file
# mesh = Mesh()
# XDMF = XDMFFile(MPI.comm_world, "traction_3Dbar.xdmf")
# XDMF.read(mesh)

# Characteristic element length - used for stabilization
h = CellDiameter(mesh)

mesh.init()
ndim = mesh.geometry().dim()  # obtain number of space dimensions
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -H/2, hsize) #and between(x[0], (0.0, 2.5))

class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H/2, hsize)

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        hsize2 = 0.02
        if near(x[0], 0.0, hsize2) and near(x[1], -H/2, hsize2) and near(x[2], T, 0.03):
            return True
        # if near(x[0], 0.0, 2*hsize) and near(x[1], H/2, 2*hsize) and near(x[2], T/2, hsize):
        #     return True

# Convert all boundary classes for visualization
bot_boundary = bot_boundary()
top_boundary = top_boundary()
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

# Constitutive functions of the damage model
def w(alpha):
    return alpha

def a(alpha):
    return (1.0-alpha)**2

def b(alpha):
    return (1.0-alpha)**6

def b_sq(alpha):
    return (1.0-alpha)**3

# Space and function, trial, and test definitions
# --------------------------------------------------------------------
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
# Originally Taylor-Hood, stabilization drops order of VectorFunctionSpace
P2 = VectorFunctionSpace(mesh, "Lagrange", 1)
P1 = FunctionSpace(mesh, "Lagrange", 1)
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
TH  = MixedElement([P2elem,P1elem])
# Mixed function space for elasticity problem (displacement and pressure) in V_u
V_u     = FunctionSpace(mesh, TH)
# Function space for damage in V_alpha
V_alpha = FunctionSpace(mesh, "Lagrange", 1)

# Define the function, trial, and test fields for elasticity problem
w_p    = Function(V_u)
u_p    = TrialFunction(V_u)
v_q    = TestFunction(V_u)
(u, p) = split(w_p)         # Displacement and pressure (u,p)
(v, q) = split(v_q)         # Test functions for (u,p)
# Define the function, test and trial fields for damage problem
alpha   = Function(V_alpha)
dalpha  = TrialFunction(V_alpha)
beta    = TestFunction(V_alpha)

# Dirichlet boundary condition (BC)
# --------------------------------------------------------------------
u00 = Expression("0.0", degree=0)
u0 = Expression(["0.0","0.0","0.0"], degree=0)
u1 = Expression("t", t=0.0, degree=0)
u2 = Expression("-t", t=0.0, degree=0)

# Displacement pinned at one point
bc_u0 = DirichletBC(V_u.sub(0).sub(2), u00, pin_point, method='pointwise')
# Displacement fixed in x and y on bottom and top boundaries
bc_u1 = DirichletBC(V_u.sub(0).sub(0), u00, bot_boundary)
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u00, bot_boundary)
bc_u3 = DirichletBC(V_u.sub(0).sub(0), u00, top_boundary)
bc_u4 = DirichletBC(V_u.sub(0).sub(1), u1, top_boundary)
bc_u = [bc_u0, bc_u1, bc_u2, bc_u3, bc_u4]

# Damage BCs to prevent damage initiating from edges
bc_alpha0 = DirichletBC(V_alpha, 0.0, bot_boundary)
bc_alpha1 = DirichletBC(V_alpha, 0.0, top_boundary)
bc_alpha = [bc_alpha0, bc_alpha1]

# Kinematics
# --------------------------------------------------------------------
d = len(u)                  # Displacement
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green (CG) tensor

# Invariants
Ic = tr(C)                      # Invariant 1 of the right CG
I2 = ((tr(C))**2-tr(C*C))/2.    # Invariant 2 of the right CG
J  = det(F)                     # Invariant 3 of the deformation gradient, F
# I3 = det(C)=J**2

# Define some parameters for the eigenvalues
d_par = tr(C)/3.
e_par = sqrt(tr((C-d_par*I)*(C-d_par*I))/6.)
# conditional(condition, true_value, false_value) where eq means ==
f_par = conditional(eq(tr((C-d_par*I)*(C-d_par*I)),0.), 0.*I, (1./e_par)*(C-d_par*I))

# IMPORTANT: Bound the argument of 'acos' both from above and from below
g_par0 = det(f_par)/2.
# ge(a,b) is a >= b and le(a,b) is a <= b
g_par1 = conditional(ge(g_par0,  1.-DOLFIN_EPS),  1.-DOLFIN_EPS, g_par0)
g_par  = conditional(le(g_par0, -1.+DOLFIN_EPS), -1.+DOLFIN_EPS, g_par1)
h_par = acos(g_par)/3.

# Define the eigenvalues of C (principal stretches) where lmbda1s <= lmbda2s<= lmbda3s
lmbda1s = d_par + 2.*e_par*cos(h_par + 2.*pi/3.)
lmbda2s = d_par + 2.*e_par*cos(h_par + 4.*pi/3.)
lmbda3s = d_par + 2.*e_par*cos(h_par + 6.*pi/3.)

# Define the energy functional of elastic problem
# --------------------------------------------------------------------
# Nominal stress tensor
def P(u, alpha):
    return a(alpha)*mu*(F - inv(F.T)) - b_sq(alpha)*p*J*inv(F.T)

# Define the Heaviside step function: H(x) = 1, x>= 0; H(x) = 0, x< 0;
def Heaviside(x):
    return conditional(ge(x, 0.), 1., 0.)

# Elastic energy, additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy  = (a(alpha)+k_ell)*(mu/2.0)*(Ic-3.0-2.*ln(J))*dx \
                - b_sq(alpha)*p*(J-1)*dx - 1./(2.*kappa)*p**2*dx
elastic_energy2 = (a(alpha)+k_ell)*(mu/2.0)*(lmbda1s+lmbda2s+lmbda3s-3.-2.*ln(lmbda1s*lmbda2s*lmbda3s))*dx \
                  - b(alpha)*p*(J-1)*dx-1./(2.*kappa)*p**2*dx

body_force        = Constant((0., 0., 0.))
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Elastic energy decomposition into active
W_act = (a(alpha)+k_ell)*Heaviside(sqrt(lmbda1s)-1.)*(mu/2.)*(lmbda1s-1.-2*ln(sqrt(lmbda1s)))*dx \
      + (a(alpha)+k_ell)*Heaviside(sqrt(lmbda2s)-1.)*(mu/2.)*(lmbda2s-1.-2*ln(sqrt(lmbda2s)))*dx \
      + (a(alpha)+k_ell)*Heaviside(sqrt(lmbda3s)-1.)*(mu/2.)*(lmbda3s-1.-2*ln(sqrt(lmbda3s)))*dx \
      + Heaviside(J-1.)*(b_sq(alpha)*p*(J-1.)+p**2/(2*kappa))
      # + b(alpha)*Heaviside(J-1.)*(kappa/2)*(J-1.)**2*dx

# Elastic energy decomposition into passive
W_pas = (mu/2.)*Heaviside(1.-sqrt(lmbda1s))*(lmbda1s-1.-ln(sqrt(lmbda1s)))*dx \
      + (mu/2.)*Heaviside(1.-sqrt(lmbda2s))*(lmbda2s-1.-ln(sqrt(lmbda2s)))*dx \
      + (mu/2.)*Heaviside(1.-sqrt(lmbda3s))*(lmbda3s-1.-ln(sqrt(lmbda3s)))*dx \
      + Heaviside(1.-J)*(p*(J-1.)+p**2/(2*kappa))
      # + Heaviside(1.-J)*(kappa/2)*(J-1.)**2*dx

# Non-dimension non-negative stability parameter
varpi_ = 1.0
# Stabilization constant
varpi = project(varpi_*h**2/(2.0*mu), FunctionSpace(mesh,'DG',0))
# Compute directional derivative about w_p in the direction of v (Gradient)
F_u = derivative(elastic_potential, w_p, v_q) \
    - varpi*b_sq(alpha)*J*inner(inv(C),outer(b_sq(alpha)*grad(p),grad(q)))*dx
    - varpi*b_sq(alpha)*J*inner(inv(C),outer(grad(b_sq(alpha))*p,grad(q)))*dx \
    + varpi*b_sq(alpha)*inner(mu*(F-inv(F.T))*grad(a(alpha)),inv(F.T)*grad(q))*dx
    # - varpi*b_sq(alpha)*J*inner(inv(C), outer(grad(p),grad(q)))*dx

# Compute directional derivative about w_p in the direction of u_p (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem for the displacement
problem_u = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solvers
solver_u  = NonlinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
# info(solver_u.parameters, True)

# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Initializing known alpha
alpha_0 = interpolate(Expression("0.", degree=0), V_alpha)
# Define the specific energy dissipation per unit volume
z = sympy.Symbol("z", positive=True)
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*dot(grad(alpha), grad(alpha)))*dx
# The elastic potential now consists of W_act + W_pas
damage_functional = W_act + W_pas + dissipated_energy

# Compute directional derivative about alpha in the direction of beta (Gradient)
E_alpha       = derivative(damage_functional, alpha, beta)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
E_alpha_alpha = derivative(E_alpha, alpha, dalpha)

# Set the lower and upper bound of the damage variable (0-1)
# initial (known) alpha lower bound
# alpha_lb = interpolate(Expression("x[0]>=-L/2 & x[0]<=0.0 & near(x[1], 0.0, 0.1*hsize) ? 1.0 : 0.0", \
#                        hsize = hsize, L=L, degree=0), V_alpha)
alpha_lb = interpolate(Expression("0.", degree=0), V_alpha)
alpha_ub = interpolate(Expression("1.", degree=0), V_alpha)  # upper bound, set to 1

# Variational problem for the damage
# (non-linear - use variational inequality solvers of petsc)
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
# info(solver_alpha.parameters,True) # uncomment to see available parameters

# Loading and initialization of vectors to store time data
load_multipliers  = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies          = np.zeros((len(load_multipliers), 5))
iterations        = np.zeros((len(load_multipliers), 2))
# Exponential function loading
# load_multipliers  = np.linspace(userpar["load_min"], userpar["load_steps"], userpar["load_steps"])
# fcn_load = []
# for steps in load_multipliers:
#     exp1 = 2.824*exp(0.001431*steps) - 2.824*exp(-0.262*steps)
#     # exp1 = 0.3304*exp(0.00189*steps) - 0.3304*exp(-0.1505*steps)
#     fcn_load.append(exp1)

# Data file name
file_tot = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
# Saves the file in case of interruption
file_tot.parameters["rewrite_function_mesh"] = False
file_tot.parameters["functions_share_mesh"]  = True
file_tot.parameters["flush_output"]          = True
# write the parameters to file
File(savedir+"/parameters.xml") << userpar

timer0 = time.process_time()        # Timer start

# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    # Printouts
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))

    # Update the displacement with each iteration
    u1.t = t

    # Alternate Mininimization scheme
    # -------------------------------------------------------------------------
    # Solve for u holding alpha constant then solve for alpha holding u constant
    iteration = 1           # Initialization of iteration loop
    err_alpha = 1.0         # Initialization for condition for iteration

    # Conditions for iterations
    while err_alpha > AM_tolerance and iteration < maxiteration:
        # solve elastic problem
        solver_u.solve()
        # solve damage problem with constraints
        solver_alpha.solve(DamageProblem(), alpha.vector(), alpha_lb.vector(), alpha_ub.vector())
        # Update the alpha condition for iteration by calculating the alpha error norm
        alpha_error = alpha.vector() - alpha_0.vector()
        err_alpha = alpha_error.norm('linf')
        # Printouts to monitor the results and number of iterations
        if MPI.rank(MPI.comm_world) == 0:
          print ("AM Iteration: {0:3d},  alpha_error: {1:>14.8f}".format(iteration, err_alpha))
        # Update variables for next iteration
        alpha_0.assign(alpha)
        iteration = iteration + 1
    # Updating the lower bound to account for the irreversibility of damage
    alpha_lb.vector()[:] = alpha.vector()

    # Split into displacement and pressure
    (u, p)   = w_p.split()

    # Project nominal stress to tensor function space
    PTensor = project(P(u, alpha), TT)
    JScalar = project(J, P1)

    # Rename for output in Paraview
    alpha.rename("Damage", "alpha")
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    PTensor.rename("Nominal Stress", "P")
    JScalar.rename("J", "J")

    # Write solution to file
    file_tot.write(alpha, t)
    file_tot.write(u, t)
    file_tot.write(p, t)
    file_tot.write(PTensor,t)
    file_tot.write(JScalar,t)

    # Post-processing
    # ----------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value  = assemble(elastic_energy)
    elastic_energy_value2 = assemble(elastic_energy2)
    surface_energy_value  = assemble(dissipated_energy)
    volume_ratio          = assemble(J/(T*H*W)*dx)

    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, elastic_energy_value+\
    	                      surface_energy_value, volume_ratio])

    # Final print outs and saving values into text files
    if MPI.rank(MPI.comm_world) == 0:
        print("\nEnd of timestep {0:3d} with load multiplier {1:4f}".format(i_t, t))
        print("\nElastic and Surface Energies: [{0:6f},{1:6f}]".format(elastic_energy_value, surface_energy_value))
        # print("\nElastic and Surface Energies: [{},{}]".format(elastic_energy_value, surface_energy_value))
        print("-----------------------------------------")
        print("\nElastic Energies [1,2]: [{},{}]".format(elastic_energy_value, elastic_energy_value2))
        print("\nVolume Ratio: [{}]".format(volume_ratio))
        # print("-----------------------------------------")
        # Save some global quantities as a function of the time
        np.savetxt(savedir + '/Taylor-Hood-energies.txt', energies)
        np.savetxt(savedir + '/Taylor-Hood-iterations.txt', iterations)

print("elapsed CPU time: ", (time.process_time() - timer0))

# Plot energy and stresses
# ----------------------------------------------------------------------------
if MPI.rank(MPI.comm_world) == 0:
    p1, = plt.plot(energies[slice(None), 0], energies[slice(None), 1])
    p2, = plt.plot(energies[slice(None), 0], energies[slice(None), 2])
    p3, = plt.plot(energies[slice(None), 0], energies[slice(None), 3])
    plt.legend([p1, p2, p3], ["Elastic", "Dissipated", "Total"], loc="best", frameon=False)
    plt.xlabel('Displacement')
    plt.ylabel('Energies')
    plt.title('Taylor-Hood FEM')
    plt.savefig(savedir + '/Taylor-Hood-energies.pdf', transparent=True)
    plt.close()
    p4, = plt.plot(energies[slice(None), 0], energies[slice(None), 4])
    plt.xlabel('Displacement')
    plt.ylabel('Volume ratio')
    plt.title('Taylor-Hood FEM')
    plt.savefig(savedir + '/Taylor-Hood-volume-ratio.pdf', transparent=True)
    plt.close()
