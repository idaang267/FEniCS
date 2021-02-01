# FEniCS code  Fracture Mechanics
################################################################################
#
# A mixed finite element method for brittle fracture in incompressible
# hyperelastic materials
#
# Modified for plane stress - no phase field terms
#
# Author: Bin Li
# Email: bl736@cornell.edu
#
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
from ufl import rank

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
                                         "maximum_iterations": 100,
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

# Initial condition (IC) class
class InitialConditions(UserExpression):
    def eval(self, values, x):
        # Displacement u0 = (values[0], values[1])
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0        # Pressure
        values[3] = 1.0        # F_33
    def value_shape(self):
         return (4,)

# Set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("mu", 1)           # Shear modulus - normalized by n*k_b*T ?
userpar.add("nu", 0.49995)     # Poisson's Ratio for slight compressibility
userpar.add("Gc", 2.4E6)       # Fracture toughness (2.4E3)
userpar.add("k_ell", 5.e-5)    # Residual stiffness
userpar.add("load_min", 0.)
userpar.add("load_max", 0.01)
userpar.add("load_steps", 5)
userpar.add("a_par", 0.0225)
userpar.add("c_par", 1.55)

# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
load_max = userpar["load_max"]
load_steps = userpar["load_steps"]

# Geometry parameters
R = 0.1                         # Radius
a_par = userpar["a_par"]
c_par = userpar["c_par"]
S = userpar["load_steps"]

# Material model parameters for plane stress
mu    = float(userpar["mu"])    # Shear modulus
nu    = userpar["nu"]           # Poisson's Ratio
E     = 2.0*(1.0+nu)*mu         # Young's Modulus
lmbda = E*nu/((1.0-nu)**2)      # Lame Parameter
kappa = (3.0-nu)/(1.0+nu)       # Bulk Modulus

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]

# Naming parameters for saving output
modelname = "2D-test"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "R_%.1f_s_%.0f_D_%.2f" % (R, load_steps, load_max)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

# For parallel processing - write one directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation of structured and refined mesh
mesh = Mesh("HalfCircleCrack.xml")
# mesh = RectangleMesh(Point(-L/2, 0), Point(L/2, H), Nx, Ny)
# Mesh printout
geo_mesh = XDMFFile(MPI.comm_world, savedir + meshname)
geo_mesh.write(mesh)

# Obtain number of space dimensions
mesh.init()
ndim = mesh.geometry().dim()
# Structure used for one printout of the statement
if MPI.rank(MPI.comm_world) == 0:
    print ("Mesh Dimension: {0:2d}".format(ndim))

# Numerical parameters of the alternate minimization scheme
maxiteration = 2000         # Sets a limit on number of iterations
AM_tolerance = 1e-4

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
            if near(x[1], 0.0, 0.001) and between(x[0],(-R+0.005, R-0.005)):
                return False
            else:
                return True
        # return on_boundary and near(x[1], 1.0, 0.01)

class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, 0.001) and between(x[0], (0.0, R))

class bot_left_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, 0.001) and between(x[0], (-R, 0.0))

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, 1E-6) and near(x[1], 0.0, 1E-6)

# Convert all boundary classes for visualization
top_boundary = top_boundary()
bot_boundary = bot_boundary()

pin_point = pin_point()
# Define lines and points
lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

# show lines of interest
lines.set_all(0)
top_boundary.mark(lines, 1)
bot_boundary.mark(lines, 2)
file_results = XDMFFile(savedir + "/" + "lines.xdmf")
file_results.write(lines)

# Show points of interest
points.set_all(0)
pin_point.mark(points, 1)
file_results = XDMFFile(savedir + "/" + "points.xdmf")
file_results.write(points)

# Variational formulation
# ----------------------------------------------------------------------------
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
# Taylor-Hood space for incompressible elasticity
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
P1 = FunctionSpace(mesh, "Lagrange", 1)
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
TH  = MixedElement([P2elem,P1elem,P1elem])
# Define function spaces for displacement, pressure, and F_{33} in V_u
V_u = FunctionSpace(mesh, TH)
# Define function space for damage in V_alpha
V_alpha = FunctionSpace(mesh, "Lagrange", 1)

# Define the function, test and trial fields for elasticity problem
w_p = Function(V_u)
u_p = TrialFunction(V_u)
v_q = TestFunction(V_u)
(u, p, F33) = split(w_p)     # Displacement, pressure, (u, p, F_{33})
(v, q, v_F33) = split(v_q)   # Test functions for u, p and F33

# Dirichlet boundary condition
# --------------------------------------------------------------------
u00 = Constant((0.0))
u0 = Expression(["0.0", "0.0"], degree=0)

# Cartesian to Polar
r_bc = Expression("sqrt(x[0]*x[0]+x[1]*x[1])",degree=1)
th_bc = Expression("atan2(x[1], x[0])", degree=1)
# Displacement x and y
u_bc = Expression(["c_par*r_bc*cos(th_bc)",
                   "a_par*sqrt(r_bc)*sin(th_bc/2)"],
                   a_par=a_par, c_par=c_par, r_bc=r_bc, th_bc=th_bc, degree=1)

# top boundary is subjected to a displacement in the y direction
bc_u1 = DirichletBC(V_u.sub(0), u_bc, top_boundary)
# bottom boundary has rollers
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u00, bot_boundary)
# Pinned center
bc_u0 = DirichletBC(V_u.sub(0), u0, pin_point, method="pointwise")

# Combined
bc_u = [bc_u0, bc_u1, bc_u2]

# Initial Conditions (IC)
#------------------------------------------------------------------------------
# Initial conditions are created by using the class defined and then
# interpolating into a finite element space
init = InitialConditions(degree=1)          # Expression requires degree def.
w_p.interpolate(init)                         # Interpolate current solution

# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
J = det(F)*F33
Ic = tr(C) + F33**2

# Define the energy functional of the elasticity problem
# --------------------------------------------------------------------
# 1st PK stress
def P(u):
    return mu*(F - inv(F.T)) - p*J*inv(F.T)

# Zero body force
body_force = Constant((0., 0.))
# Elastic energy, additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy    = (mu/2.)*(Ic-3.-2.*ln(J))*dx \
                    - p*(J-1.)*dx - 1./(2.*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Compute directional derivative about w_p in the direction of v (Gradient)
# Define the additional weak form eq. for plane stress
F_u = derivative(elastic_potential, w_p, v_q) \
     + (F33**2 - 1. - p*J/mu)*v_F33*dx
# Compute directional derivative about w_p in the direction of u_p (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem to solve for displacement and pressure
problem_up = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solver for displacement and pressure
solver_up  = NonlinearVariationalSolver(problem_up)
solver_up.parameters.update(solver_up_parameters)
# info(solver_up.parameters, True) # uncomment to see available parameters

# loading and initialization of vectors to store time datas
a_par_list = np.linspace(userpar["load_min"], userpar["a_par"], userpar["load_steps"])
c_par_list = np.linspace(userpar["load_min"], userpar["c_par"], userpar["load_steps"])
# energies       = np.zeros((len(load_multipliers), 5))
# iterations     = np.zeros((len(load_multipliers), 2))

# Split into displacement and pressure
(u, p, F33) = w_p.split()
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
for (i_t, t) in enumerate(a_par_list):
    a_par = a_par_list[i_t]
    c_par = c_par_list[i_t]

    # Update the displacement with each iteration
    u_bc.a_par = a_par
    u_bc.c_par = c_par

    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of step {0:2d}: a_par = {1:4f} c_par = {2:4f} ---\033[1;m".format(i_t, a_par, c_par))

    # Solve elastic problem
    solver_up.solve()

    # Project nominal stress to tensor function space
    PTensor = project(P(u), TT)
    JScalar = project(J, P1)

    # Rename for paraview
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    F33.rename("F33", "F33")
    PTensor.rename("Nominal Stress", "P")
    JScalar.rename("J", "J")

    # Write solution to file
    file_tot.write(u, t)
    file_tot.write(p, t)
    file_tot.write(F33, t)
    file_tot.write(PTensor,t)
    file_tot.write(JScalar,t)

# ----------------------------------------------------------------------------
print("elapsed CPU time: ", (time.process_time() - timer0))

'''
# Plot energy and stresses
if MPI.rank(MPI.comm_world) == 0:
    p1, = plt.plot(energies[slice(None), 0], energies[slice(None), 1])
    p2, = plt.plot(energies[slice(None), 0], energies[slice(None), 2])
    p3, = plt.plot(energies[slice(None), 0], energies[slice(None), 3])
    plt.legend([p1, p2, p3], ["Elastic", "Dissipated", "Total"], loc="best", frameon=False)
    plt.xlabel('Displacement')
    plt.ylabel('Energies')
    plt.title('stabilized FEM')
    plt.savefig(savedir + '/stabilized-energies.pdf', transparent=True)
    plt.close()
    p4, = plt.plot(energies[slice(None), 0], energies[slice(None), 4])
    plt.xlabel('Displacement')
    plt.ylabel('Volume ratio')
    plt.title('stabilized FEM')
    plt.savefig(savedir + '/stabilized-volume-ratio.pdf', transparent=True)
    plt.close()
'''
