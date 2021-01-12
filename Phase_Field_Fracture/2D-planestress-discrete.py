# FEniCS code  Fracture Mechanics
################################################################################
#
# A stabilized mixed finite element method for brittle fracture in
# incompressible hyperelastic materials
#
# Modified for plane stress cases
#
# Author: Bin Li
# Email: bl736@cornell.edu
# date: 10/01/2018
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
userpar.add("meshsizeX", 100)
userpar.add("meshsizeY", 100)
userpar.add("load_min", 0.)
userpar.add("load_max", 1.0)
userpar.add("load_steps", 10)
userpar.add("KI",1.0)        # mode I loading

# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
load_max = userpar["load_max"]
load_steps = userpar["load_steps"]

# Geometry parameters
L, H = 1.0, 1.0        # Length (x) and height (y-direction)
Nx   = userpar["meshsizeX"]
Ny   = userpar["meshsizeY"]
hsize = float(H/Ny)    # Geometry based definition for regularization
S = userpar["load_steps"]

# Material model parameters
mu    = float(userpar["mu"])           # Shear modulus
nu    = userpar["nu"]           # Poisson's Ratio
lmbda = 2.0*mu*nu/(1.0-2.0*nu)  # Lame Parameter
kappa = 2*(1+nu)*mu/(3*(1-2*nu))# Bulk Modulus

# Fracture toughness and residual stiffness
KI    = float(userpar["KI"])
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]

# Naming parameters for saving output
modelname = "2D-test"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "Nx_%.0f_Ny_%.0f_s_%.0f" % (Nx, Ny, load_steps)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

# For parallel processing - write one directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation of structured and refined mesh
mesh = Mesh("SimpleCrack.xml")
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

# Reference value for the loading (imposed displacement)
ut = 1.0

# Numerical parameters of the alternate minimization scheme
maxiteration = 2000         # Sets a limit on number of iterations
AM_tolerance = 1e-4

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class right_top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.5, 0.1*hsize) and between(x[1], (0.0, 0.5))

class right_bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.5, 0.1*hsize) and between(x[1], (-0.5, 0.0))

class left_top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -0.5, 0.1*hsize) and between(x[1], (0.0, 0.5))

class left_bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -0.5, 0.1*hsize) and between(x[1], (-0.5, 0.0))

class bot_left_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -0.5, 0.1 *hsize) and between(x[0], (-0.5, -0.01))

class bot_right_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], -0.5, 0.1*hsize) and between(x[0], (0.0, 0.5))

class top_left_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.5, 0.1*hsize) and between(x[0], (-0.5, -0.005))

class top_right_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.5, 0.1*hsize) and between(x[0], (0, 0.5))

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, 0.01) and near(x[1], 0.0, 0.01)

# Convert all boundary classes for visualization
right_top_boundary = right_top_boundary()
right_bot_boundary = right_bot_boundary()
left_top_boundary = left_top_boundary()
left_bot_boundary = left_bot_boundary()
bot_left_boundary = bot_left_boundary()
bot_right_boundary = bot_right_boundary()
top_left_boundary = top_left_boundary()
top_right_boundary = top_right_boundary()
pin_point = pin_point()
# Define lines and points
lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

# show lines of interest
lines.set_all(0)
right_top_boundary.mark(lines, 1)
right_bot_boundary.mark(lines, 2)
left_top_boundary.mark(lines, 3)
left_bot_boundary.mark(lines, 4)
bot_left_boundary.mark(lines, 5)
bot_right_boundary.mark(lines, 6)
top_left_boundary.mark(lines, 7)
top_right_boundary.mark(lines, 8)
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
# dim = V_u.dim()
# ndim = mesh.geometry().dim()
# coords = V_u.tabulate_dof_coordinates().reshape(dim, ndim)
# xcoords = coords[:,0]
# ycoords = coords[:,1]
# POI_1 = np.where(xcoords == 0.5)[0]
# POI_2 = np.where(np.logical_and(ycoords >= -0.5, ycoords <= 0.5))[0]
# cor_RHS = coords[POI_1]
# np.savetxt(savedir + '/coordinates.txt', cor_RHS)

u00 = Constant((0.0))
u0 = Expression(["0.0", "0.0"], degree=0)

# Cartesian to Polar
r_bc = Expression("sqrt(x[0]*x[0]+x[1]*x[1])",degree=1)
th_r = Expression("atan(x[1]/x[0])", degree=1)
th_t = Expression("atan(x[1]/x[0]) + pi", degree=1)
th_b = Expression("atan(x[1]/x[0]) - pi", degree=1)

ramp_par = 0.0
# Displacement x and y
u_r = Expression(["ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(7/3*cos(th_r/2) - cos(3*th_r/2))) + \
                   ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(-1/15*std::log(r_bc)-52/45*(std::log(r_bc)+3/4*sin(th_r)*sin(th_r)) - 103/48*cos(th_r) + 26/15*cos(2*th_r)-3/16*cos(3*th_r))",
                  "ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(13/3*sin(th_r/2) - sin(3*th_r/2))) + \
                  ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(th_r/15 - 52/45*(th_r/4 - 3/8*sin(2*th_r)) - 61/48*sin(th_r) + 26/15*sin(2*th_r) - 3/16*sin(3*th_r))"],
                  ramp_par=ramp_par, KI=KI, mu=mu, r_bc=r_bc, th_r=th_r, degree=1)

u_t = Expression(["ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(7/3*cos(th_t/2) - cos(3*th_t/2))) + \
                   ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(-1/15*std::log(r_bc)-52/45*(std::log(r_bc)+3/4*sin(th_t)*sin(th_t)) - 103/48*cos(th_t) + 26/15*cos(2*th_t)-3/16*cos(3*th_t))",
                  "ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(13/3*sin(th_t/2) - sin(3*th_t/2))) + \
                   ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(th_t/15 - 52/45*(th_t/4 - 3/8*sin(2*th_t)) - 61/48*sin(th_t) + 26/15*sin(2*th_t) - 3/16*sin(3*th_t))"],
                  ramp_par=ramp_par, KI=KI, mu=mu, r_bc=r_bc, th_t=th_t, degree=1)

u_b = Expression(["ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(7/3*cos(th_b/2) - cos(3*th_b/2))) + \
                   ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(-1/15*std::log(r_bc)-52/45*(std::log(r_bc)+3/4*sin(th_b)*sin(th_b)) - 103/48*cos(th_b) + 26/15*cos(2*th_b)-3/16*cos(3*th_b))",
                  "ramp_par*(KI/(4*mu))*(sqrt(r_bc/(2*pi))*(13/3*sin(th_b/2) - sin(3*th_b/2))) + \
                   ramp_par*ramp_par*(KI/(4*mu*sqrt(2*pi)))*(KI/(4*mu*sqrt(2*pi)))*(th_b/15 - 52/45*(th_b/4 - 3/8*sin(2*th_b)) - 61/48*sin(th_b) + 26/15*sin(2*th_b) - 3/16*sin(3*th_b))"],
                  ramp_par=ramp_par, KI=KI, mu=mu, r_bc=r_bc, th_b=th_b, degree=1)

# Top boundary
bc_u3 = DirichletBC(V_u.sub(0), u_r, top_right_boundary)
bc_u4 = DirichletBC(V_u.sub(0), u_t, top_left_boundary)
# # Right boundary
bc_u1 = DirichletBC(V_u.sub(0), u_r, right_top_boundary)
bc_u2 = DirichletBC(V_u.sub(0), u_r, right_bot_boundary)
# # Bottom boundary
bc_u5 = DirichletBC(V_u.sub(0), u_r, bot_right_boundary)
bc_u6 = DirichletBC(V_u.sub(0), u_b, bot_left_boundary)
# Left boundary
bc_u7 = DirichletBC(V_u.sub(0), u_t, left_top_boundary)
bc_u8 = DirichletBC(V_u.sub(0), u_b, left_bot_boundary)

bc_u9 = DirichletBC(V_u.sub(0), u0, pin_point)

# Combine
bc_u = [bc_u1, bc_u2, bc_u3, bc_u4, bc_u5, bc_u6, bc_u7, bc_u8]

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
elastic_energy    = (mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
                    - p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Define the stabilization term and the additional weak form eq.
# Compute directional derivative about w_p in the direction of v (Gradient)
F_u = derivative(elastic_potential, w_p, v_q) \
     + (F33**2 - 1 + p*J/mu)*v_F33*dx
# Compute directional derivative about w_p in the direction of u_p (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem to solve for displacement and pressure
problem_up = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solver for displacement and pressure
solver_up  = NonlinearVariationalSolver(problem_up)
solver_up.parameters.update(solver_up_parameters)
# info(solver_up.parameters, True) # uncomment to see available parameters

# loading and initialization of vectors to store time datas
test_scale = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
# energies         = np.zeros((len(load_multipliers), 5))
# iterations       = np.zeros((len(load_multipliers), 2))

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
for (i_t, t) in enumerate(test_scale):
    ramp_par = test_scale[i_t]

    # Update the displacement with each iteration
    u_r.ramp_par = ramp_par
    u_t.ramp_par = ramp_par
    u_b.ramp_par = ramp_par

    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of step {0:2d}: test = {1:2f} ---\033[1;m".format(i_t, test_scale[i_t]))

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
