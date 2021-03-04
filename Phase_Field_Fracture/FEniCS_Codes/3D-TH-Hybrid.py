# -------------------------------------------
# FEniCS code  Variational Fracture Mechanics
################################################################################
#                                                                              #
# A Taylor-Hood finite element method for gradient damage models of            #
# fracture in incompressible hyperelastic materials                            #
# author: Bin Li                                                               #
# Email: bl736@cornell.edu                                                     #
# date: 10/01/2018                                                             #
#                                                                              #
################################################################################
# e.g. python3 traction-neo-Hookean.py --meshsize 100						   #
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
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# parameters of the solvers
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

# parameters of the PETSc/Tao solver used for the alpha-problem
tao_solver_parameters = {"maximum_iterations": 100,
                         "report": False,
                         "line_search": "more-thuente",
                         "linear_solver": "cg",
                         "preconditioner" : "hypre_amg",
                         "method": "tron",
                         "gradient_absolute_tol": 1e-8,
                         "gradient_relative_tol": 1e-8,
                         "error_on_nonconvergence": True}

# Implement the box constraints for damage field
# --------------------------------------------------------------------
# Variational problem for the damage (non-linear to use variational inequality solvers of petsc)
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

# User Parameters
# ----------------------------------------------------------------------------
parameters.parse()
userpar = Parameters("user")
userpar.add("mu", 1)          # Shear modulus
userpar.add("kappa",1000)
userpar.add("Gc",2.4e6)       # fracture toughness
userpar.add("k_ell",5.e-5)    # residual stiffness
userpar.add("alpha_t", 0.8)   # Damage threshold
userpar.add("meshsizeX", 40)
userpar.add("meshsizeY", 18)
userpar.add("meshsizeZ", 18)
userpar.add("load_min",0.)
userpar.add("load_max", 0.01)
userpar.add("load_steps", 4)
userpar.add("ell_multi", 2)
# Parse command-line options
userpar.parse()

# Geometry parameters
L, H, W = 0.5, 0.1, 0.1
Nx   = userpar["meshsizeX"]
Ny   = userpar["meshsizeY"]
Nz   = userpar["meshsizeZ"]
# Geometry based definition of element size in the y direction
hsize   = float(H/Ny)

# Material constants
mu    = userpar["mu"]
kappa = userpar["kappa"]
alpha_t = userpar["alpha_t"]

# Fracture toughness and residual stiffness
Gc    = userpar["Gc"]
k_ell = userpar["k_ell"]
ell_multi = userpar["ell_multi"]
# regularization paramater
ell = Constant(ell_multi*hsize)

# Numerical parameters of the alternate minimization
maxiteration = 2000
AM_tolerance = 1e-4

modelname = "3D-TaylorHood"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "Hybrid_Nx_%.0f_Ny_%.0f_Nz_%.0f_ellx_%.0f_at_%.2f" % (Nx, Ny, Nz, ell_multi, alpha_t)
savedir   = "output/" + modelname + "/" + simulation_params + "/"

if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(L, H, W), Nx, Ny, Nz)
# mesh = Mesh("ShearTestRefStruct.xml")
geo_mesh  = XDMFFile(MPI.comm_world, savedir+meshname)
geo_mesh.write(mesh)

mesh.init()
ndim = mesh.geometry().dim()  # get number of space dimensions
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
class bot_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0, hsize)

class top_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], H, hsize)

class back_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 0, hsize)

class front_boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], W, hsize)

class pin_point(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L/2, hsize) and near(x[1], H/2, hsize) and near(x[2], W/2, hsize)

# Convert all boundary classes for visualization
bot_boundary = bot_boundary()
top_boundary = top_boundary()
back_boundary = back_boundary()
front_boundary = front_boundary()
pin_point = pin_point()

lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 3)

# show lines of interest
lines.set_all(0)
bot_boundary.mark(lines, 1)
top_boundary.mark(lines, 1)
file_results = XDMFFile(savedir + "/" + "lines.xdmf")
file_results.write(lines)

# Variational formulation
# ----------------------------------------------------------------------------
# Tensor space for projection of stress
TT = TensorFunctionSpace(mesh,'DG',0)
# Create Taylor-Hood function space for elasticity + Damage
P2 = VectorFunctionSpace(mesh, "Lagrange", 2)
P1 = FunctionSpace(mesh, "Lagrange", 1)
P2elem = P2.ufl_element()
P1elem = P1.ufl_element()
TH  = MixedElement([P2elem,P1elem])

# Define function spaces for displacement, pressure, and F_{33} in V_u
V_u     = FunctionSpace(mesh, TH)
V_alpha = FunctionSpace(mesh, "Lagrange", 1)

# Define the function, test and trial fields for elasticity problem
w_p     = Function(V_u)
u_p     = TrialFunction(V_u)
v_q     = TestFunction(V_u)
(u, p)  = split(w_p)

# define the function, test, and trial for damage problem
alpha   = Function(V_alpha)
dalpha  = TrialFunction(V_alpha)
beta    = TestFunction(V_alpha)

# Dirichlet boundary condition
# --------------------------------------------------------------------
u00 = Constant((0.0))
u0 = Expression(["0.0","0.0","0.0"], degree=0)
u1 = Expression("t",  t=0.0, degree=0)
u2 = Expression("-t", t=0.0, degree=0)
# Displacement on top and bottom boundaries
bc_u1 = DirichletBC(V_u.sub(0).sub(1), u1, top_boundary)
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u2, bot_boundary)
# Rollers on front and back boundary constraining displacement in the z
bc_u3 = DirichletBC(V_u.sub(0).sub(2), u00, back_boundary)
bc_u4 = DirichletBC(V_u.sub(0).sub(2), u00, front_boundary)
# Pinned center point
bc_u0 = DirichletBC(V_u.sub(0), u0, pin_point)
bc_u = [bc_u0, bc_u1, bc_u2, bc_u3, bc_u4]

# bc - alpha (zero damage)
bc_alpha0 = DirichletBC(V_alpha, 0.0, top_boundary)
bc_alpha1 = DirichletBC(V_alpha, 0.0, bot_boundary)
bc_alpha = [bc_alpha0, bc_alpha1]

# Constitutive functions of the damage model
# ----------------------------------------------------------------------------
def w(alpha):
    return alpha

def a(alpha):
    return (1.0-alpha)**2

def b(alpha):
    return (1.0-alpha)**3

# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Define damage criterions according to damage threshold, alpha_t
def d_c_d(alpha, alpha_t):  # elastic energy criterion
    return conditional(ge(alpha, alpha_t), 0, 1)

def d_c_e(alpha, alpha_t):  # damage degraded energy criterion
    return conditional(le(alpha, alpha_t), 0, 1)

# Nominal stress tensor
def P(u, alpha, alpha_t):
    return conditional(ge(alpha, alpha_t), 0, 1)*mu*(F - inv(F.T)) \
    - conditional(ge(alpha, alpha_t), 0, 1)*p*J*inv(F.T) \
    + conditional(le(alpha, alpha_t), 0, 1)*a(alpha)*mu*(F - inv(F.T)) \
    - conditional(le(alpha, alpha_t), 0, 1)*b(alpha)*p*J*inv(F.T)

# Elastic energy, additional terms enforce material incompressibility and regularizes the Lagrange Multiplier
elastic_energy = d_c_d(alpha, alpha_t)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
               - d_c_d(alpha, alpha_t)*p*(J-1.0)*dx \
               - d_c_e(alpha, alpha_t)*(a(alpha)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
               - d_c_e(alpha, alpha_t)*b(alpha)*p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
body_force     = Constant((0., 0., 0.))
external_work  = dot(body_force, u)*dx
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

# --------------------------------------------------------------------
# Define the energy functional of damage problem
# --------------------------------------------------------------------
alpha_0 = interpolate(Expression("0.", degree=0), V_alpha)  # initial (known) alpha
z = sympy.Symbol("z", positive=True)
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*dot(grad(alpha), grad(alpha)))*dx
damage_functional = elastic_potential + dissipated_energy

# Compute directional derivative about alpha in the direction of beta (Gradient)
E_alpha       = derivative(damage_functional, alpha, beta)
# Compute directional derivative about alpha in the direction of dalpha (Hessian)
E_alpha_alpha = derivative(E_alpha, alpha, dalpha)

# Set up the solvers
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
# info(solver_alpha.parameters,True) # uncomment to see available parameters

#alpha_lb = interpolate(Expression("x[0]<=0.5*L & near(x[1], 0.0, tol) ? 1.0 : 0.0", \
#                                  degree=0, L= L, tol=1.2*cra_w), V_alpha)  # initial (known) alpha
# alpha_lb = interpolate(Expression("0.", degree=0), V_alpha)  # lower bound, set to 0
alpha_lb = interpolate(Expression("x[0]>=0 & x[0]<=L/2 & near(x[1], H/2, 0.1 * hsize) ? 1.0 : 0.0", \
                       hsize = hsize, L=L, H=H, degree=0), V_alpha)
alpha_ub = interpolate(Expression("1.", degree=0), V_alpha)  # upper bound, set to 1

# loading and initialization of vectors to store time datas
load_multipliers = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies         = np.zeros((len(load_multipliers), 5))
iterations       = np.zeros((len(load_multipliers), 2))

# set the saved data file name
(u, p)      = w_p.split()
# Data file name
file_tot = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
# Saves the file in case of interruption
file_tot.parameters["rewrite_function_mesh"] = False
file_tot.parameters["functions_share_mesh"]  = True
file_tot.parameters["flush_output"]          = True
# write the parameters to file
File(savedir+"/parameters.xml") << userpar

# ----------------------------------------------------------------------------
# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    u1.t = t
    u2.t = t

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
        solver_alpha.solve(DamageProblem(), alpha.vector(), alpha_lb.vector(), alpha_ub.vector())
        # test error
        alpha_error = alpha.vector() - alpha_0.vector()
        err_alpha = alpha_error.norm('linf')
        # monitor the results
        volume_ratio = assemble(J/(L*H*W)*dx)
        if MPI.rank(MPI.comm_world) == 0:
          print ("AM Iteration: {0:3d},  alpha_error: {1:>14.8f}".format(iteration, err_alpha))
          print("\nVolume Ratio: [{}]".format(volume_ratio))
        # update iteration
        alpha_0.assign(alpha)
        iteration = iteration + 1
    # updating the lower bound to account for the irreversibility
    alpha_lb.vector()[:] = alpha.vector()

    # Project nominal stress to tensor function space
    PTensor = project(P(u, alpha, alpha_t), TT)
    JScalar = project(J, P1)

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

    # Some post-processing
    # ----------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value = assemble(elastic_energy)
    surface_energy_value = assemble(dissipated_energy)
    volume_ratio         = assemble(J/(L*H*W)*dx)

    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, elastic_energy_value+\
    	                      surface_energy_value, volume_ratio])

    if MPI.rank(MPI.comm_world) == 0:
        print("\nEnd of timestep {0:3d} with load multiplier {1:4f}".format(i_t, t))
        print("\nElastic and Surface Energies: [{0:6f},{1:6f}]".format(elastic_energy_value, surface_energy_value))
        print("\nElastic and Surface Energies: [{},{}]".format(elastic_energy_value, surface_energy_value))
        print("\nVolume Ratio: [{}]".format(volume_ratio))
        print("-----------------------------------------")
        # Save some global quantities as a function of the time
        np.savetxt(savedir + '/Taylor-Hood-energies.txt', energies)
        np.savetxt(savedir + '/Taylor-Hood-iterations.txt', iterations)
# ----------------------------------------------------------------------------

# Plot energy and stresses
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
