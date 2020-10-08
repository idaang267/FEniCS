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
from ufl import rank

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
userpar.add("mu", 1)          # Shear modulus - normalized by n*k_b*T ?
userpar.add("nu", 0.49995)     # Poisson's Ratio for slight compressibility
userpar.add("Gc", 2.4E6)       # Fracture toughness (2.4E3)
userpar.add("k_ell", 5.e-5)    # Residual stiffness
userpar.add("meshsizeX", 1000)
userpar.add("meshsizeY", 500)
userpar.add("load_min", 0.)
userpar.add("load_max", 1.0)
userpar.add("load_steps", 21)
userpar.add("ell_multi", 20)
# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
# Geometry parameters
L, H = 15, 0.75        # Length (x) and height (y-direction)
Nx   = userpar["meshsizeX"]
Ny   = userpar["meshsizeY"]
# Ny = int(N*H/L)
hsize = float(H/Ny)    # Geometry based definition for regularization
S = userpar["load_steps"]

# Material model parameters
mu    = userpar["mu"]           # Shear modulus
nu    = userpar["nu"]           # Poisson's Ratio
lmbda = 2.0*mu*nu/(1.0-2.0*nu)  # Lame Parameter
kappa = 2*(1+nu)*mu/(3*(1-2*nu))# Bulk Modulus

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]
ell_multi = userpar["ell_multi"]

# Naming parameters for saving output
modelname = "2D-discrete-stabilized"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "Nx_%.0f_Ny_%.0f_S_%.0f_ellx_%.0f" % (Nx, Ny, S, ell_multi)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

# For parallel processing - write one directory
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation of structured mesh
mesh = RectangleMesh(Point(0, 0), Point(L, H), Nx, Ny)
# Mesh generation of structured and refined mesh
# mesh = Mesh("2DShearTestRef.xml")
# Mesh rpintout
# geo_mesh = XDMFFile(MPI.comm_world, savedir + meshname)
# geo_mesh.write(mesh)

# Obtain number of space dimensions
mesh.init()
ndim = mesh.geometry().dim()
# Structure used for one printout of the statement
if MPI.rank(MPI.comm_world) == 0:
    print ("Mesh Dimension: {0:2d}".format(ndim))

# Damage regularization parameter - internal length scale used for tuning Gc
ell = Constant(ell_multi*hsize)
# Characteristic element length - used for stabilization
h = CellDiameter(mesh)

# Reference value for the loading (imposed displacement)
ut = 1.0

# Numerical parameters of the alternate minimization scheme
maxiteration = 2000         # Sets a limit on number of iterations
AM_tolerance = 1e-4

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class right_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, hsize)

class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and between(x[0], (5, L)) and near(x[1], 0, hsize)

class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, hsize)

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, hsize) and near(x[1], 0, hsize)

# Convert all boundary classes for visualization
right_boundary = right_boundary()
bot_boundary = bot_boundary()
top_boundary = top_boundary()
pin_point = pin_point()

lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

# show lines of interest
lines.set_all(0)
right_boundary.mark(lines,1)
bot_boundary.mark(lines, 2)
top_boundary.mark(lines, 3)
file_results = XDMFFile(savedir + "/" + "lines.xdmf")
file_results.write(lines)

# Show points of interest
points.set_all(0)
pin_point.mark(points, 1)
file_results = XDMFFile(savedir + "/" + "points.xdmf")
file_results.write(points)

# Variational formulation
# ----------------------------------------------------------------------------
P2 = VectorFunctionSpace(mesh, "Lagrange", 1)
P1 = FunctionSpace(mesh, "Lagrange", 1)
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
TH  = MixedElement([P2elem,P1elem,P1elem])
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
u1 = Expression("t", t= 0.0, degree=0)
# Pin Point
bc_u0 = DirichletBC(V_u.sub(0), u0, pin_point)
# Roller
bc_u1 = DirichletBC(V_u.sub(0).sub(1), u00, bot_boundary)
# Displacement
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u1, top_boundary)

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
# Zero body force
body_force = Constant((0., 0.))
# Non-dimension non-negative stability parameter
varpi_ = 1.0
# Eq 19 in Klaas
varpi  = project(varpi_*h**2/(2.0*mu), FunctionSpace(mesh,'DG',0))
# Elastic energy, additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy    = (mu/2.0)*(Ic-3.0-2.0*ln(J))*dx + p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Define the energy functional of the damage problem - the fracture energy
# dissipated_energy = Gc*dx

# Define the stabilization term and the additional weak form eq.
# Compute directional derivative about w_p in the direction of v (Gradient)
F_u = derivative(elastic_potential, w_p, v_q) \
        - varpi*J*inner(inv(C), outer(grad(p),grad(q)))*dx \
      + (F33**2 - 1 + p*J/mu)*v_F33*dx
# Compute directional derivative about w_p in the direction of u_p (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem to solve for displacement and pressure
problem_up = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solver for displacement and pressure
solver_up  = NonlinearVariationalSolver(problem_up)
solver_up.parameters.update(solver_up_parameters)
# info(solver_up.parameters, True) # uncomment to see available parameters

# Loading and initialization of vectors to store data of interest
load_multipliers = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies         = np.zeros((len(load_multipliers), 5))
iterations       = np.zeros((len(load_multipliers), 2))

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
for (i_t, t) in enumerate(load_multipliers):
    # Update the displacement with each iteration
    u1.t = t
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))

    # Solve elastic problem
    solver_up.solve()

    # Rename for paraview
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    F33.rename("F33", "F33")

    # Write solution to file
    file_tot.write(u, t)
    file_tot.write(p, t)
    file_tot.write(F33, t)

# ----------------------------------------------------------------------------
print("elapsed CPU time: ", (time.process_time() - timer0))

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
