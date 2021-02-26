# -------------------------------------------
# FEniCS code  Variational Fracture Mechanics
################################################################################
#                                                                              #
# A stabilized mixed finite element method for gradient damage models of 	   #
# fracture in incompressible hyperelastic materials                            #
# author: Bin Li                                                               #
# Email: bl736@cornell.edu                                                     #
# Modified for plane strain
#                                                                              #
################################################################################
# e.g. python3 traction-stabilized.py --meshsize 100						   #
################################################################################

# ----------------------------------------------------------------------------
from __future__ import division
from dolfin import *
from mshr import *
from scipy import optimize
from ufl import rank

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
# --------------------------------------------------------------------
# Parameters of the nonlinear SNES solver used for the displacement u-problem
solver_u_parameters   = {"nonlinear_solver": "snes",
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

# Set up the solvers
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
# info(solver_alpha.parameters,True) # uncomment to see available parameters

# Variational problem for the damage problem (non-linear - use variational inequality solvers in PETSc)
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
userpar.add("mu", 1)     # Shear modulus - normalized by n*k_b*T ?
userpar.add("kappa",10)
userpar.add("Gc", 2.4e6)       # Fracture toughness
userpar.add("k_ell", 5.e-5)    # residual stiffness
userpar.add("meshsizeX", 1000)
userpar.add("meshsizeY", 200)
userpar.add("load_min", 0.)
userpar.add("load_max", 0.1)
userpar.add("load_steps", 5)
# Parse command-line options
userpar.parse()

# Constants: some parsed from user parameters
# ----------------------------------------------------------------------------
# Geometry parameters
# Length (x) and height (y-direction)
L, H = 5, 1.0
Nx   = userpar["meshsizeX"]
Ny   = userpar["meshsizeY"]
hsize = float(H/Ny)    # Geometry based definition for regularization
S = userpar["load_steps"]

# Material constants parsed from user parameters
# ----------------------------------------------------------------------------
mu    = userpar["mu"]           # Shear modulus
kappa = userpar["kappa"]        # Bulk Modulus

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]

# Naming for the saved output
modelname = "2D-planestrain-discrete"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "Nx_%.0f_Ny_%.0f_d_%.2f_k_%.0f" % (Nx, Ny, userpar["load_max"], kappa)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

# For parallel processing
if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation of structured mesh
mesh = Mesh("2DShearTestHalf.xml")
# mesh = RectangleMesh(Point(-L/2, -H/2), Point(L/2, H/2), Nx, Ny)
# Mesh generation of structured and refined mesh
# mesh = Mesh("2DShearTestRef.xml")
# Mesh printout
# geo_mesh = XDMFFile(MPI.comm_world, savedir + meshname)
# geo_mesh.write(mesh)

# Obtain number of space dimensions
mesh.init()
ndim = mesh.geometry().dim()
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# Loading Parameters
ut = 1.0           # Reference value for the loading (imposed displacement)

# Numerical parameters of the alternate minimization
maxiteration = 2000
AM_tolerance = 1e-4

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1.0, 0.1*hsize)

class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, 0.1 *hsize) and between(x[0], (0.0, 2.5))

class right_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 2.5, 0.1*hsize)

class left_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], -2.5, 0.1*hsize)

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, 0.1*hsize) and near(x[1], 0.0, 0.1*hsize)

# Convert all boundary classes for visualization
top_boundary = top_boundary()
bot_boundary = bot_boundary()
right_boundary = right_boundary()
left_boundary = left_boundary()
pin_point = pin_point()

# Define lines and points
lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

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
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
# Taylor-Hood space for incompressible elasticity
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
P1 = FunctionSpace(mesh, "Lagrange", 1)
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
TH  = MixedElement([P2elem,P1elem])

# # Create function space for elasticity + damage
# P2  = VectorElement("Lagrange", mesh.ufl_cell(), 2)
# P1  = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
# # Stabilized mixed FEM for incompressible elasticity
# TH  = MixedElement([P2,P1])
# Define function spaces for displacement and pressure in V_u
V_u = FunctionSpace(mesh, TH)

# Define the function, test and trial fields for elasticity problem
w_p    = Function(V_u)
u_p    = TrialFunction(V_u)
v_q    = TestFunction(V_u)
(u, p) = split(w_p)     # Displacement and pressure
(v, q) = split(v_q)     # Test functions for displacement and pressure

# Dirichlet boundary condition
# --------------------------------------------------------------------
u00 = Constant((0.0))
u0 = Expression(["0.0", "0.0"], degree=0)
u1 = Expression("t", t= 0.0, degree=0)
u2 = Expression("-t", t= 0.0, degree=0)

# Top/bottom boundaries have displacement in the y direction
bc_u1 = DirichletBC(V_u.sub(0).sub(1), u1, top_boundary)
# bottom boundary has rollers
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u00, bot_boundary)
# Pointwise
bc_u0 = DirichletBC(V_u.sub(0), u0, pin_point, method="pointwise")

# combine boundary conditions
bc_u = [bc_u0, bc_u1, bc_u2]

# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
J  = det(F)
Ic = tr(C) + 1

# Nominal stress tensor
def P(u):
    return mu*(F - inv(F.T)) - p*J*inv(F.T)

# Zero body force
body_force        = Constant((0., 0.))
# Elastic energy, additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy    = (mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
                    - p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# Compute directional derivative about w_p in the direction of v (Gradient)
F_u = derivative(elastic_potential, w_p, v_q)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
J_u = derivative(F_u, w_p, u_p)

# Variational problem for the displacement
problem_u = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
# Set up the solvers
solver_u  = NonlinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
# info(solver_u.parameters, True)

# Loading and initialization of vectors to store time datas
load_multipliers  = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies          = np.zeros((len(load_multipliers), 5))

# Set the saved data file name
(u, p)     = w_p.split()
file_tot   = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
# Saves the file in case of interruption
file_tot.parameters["rewrite_function_mesh"] = False
file_tot.parameters["functions_share_mesh"]  = True
file_tot.parameters["flush_output"]          = True
# Write the parameters to file
File(savedir+"/parameters.xml") << userpar

timer0 = time.process_time()

# ----------------------------------------------------------------------------
# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    # Update the displacement with each iteration
    u1.t = t

    # Structure used for one printout of the statement
    if MPI.rank(MPI.comm_world) == 0:
        print("\033[1;32m--- Starting of Time step {0:2d}: t = {1:4f} ---\033[1;m".format(i_t, t))

    # Solve elastic problem
    solver_u.solve()

    # Project nominal stress to tensor function space
    PTensor = project(P(u), TT)
    JScalar = project(J, P1)

    # Rename for paraview
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    PTensor.rename("Nominal Stress", "P")
    JScalar.rename("J", "J")

    # Write solution to file
    file_tot.write(u, t)
    file_tot.write(p, t)
    file_tot.write(PTensor,t)
    file_tot.write(JScalar,t)

#----------------------------------------------------------------------
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
