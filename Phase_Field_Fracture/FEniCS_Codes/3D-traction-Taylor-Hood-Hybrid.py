# -------------------------------------------
# FEniCS code  Variational Fracture Mechanics
#################################################################################################################
#                                                                                                               #
# A Taylor-Hood finite element method for gradient damage models of                                             #
# fracture in incompressible hyperelastic materials                                                             #
# author: Bin Li                                                                                                #
# Email: bl736@cornell.edu                                                                                      #
# date: 10/01/2018                                                                                              #
#                                                                                                               #
#################################################################################################################
# e.g. python3 traction-neo-Hookean.py --meshsize 100														                              	#
#################################################################################################################


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
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
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

# set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("mu",5.0e4)       # n*k*T Shear modulus
userpar.add("Gc",2.4e3)       # fracture toughness
userpar.add("nu",0.4995)      # bulk modulus for slightly compressibility
userpar.add("k_ell",5.e-5)    # residual stiffness
userpar.add("meshsize",100)
userpar.add("load_min",0.00)
userpar.add("load_max",0.50)
userpar.add("load_steps",101)
# Parse command-line options
userpar.parse()

# Material constants
# ----------------------------------------------------------------------------
mu    = userpar["mu"]
Gc    = userpar["Gc"]
nu    = userpar["nu"]
k_ell = userpar["k_ell"]

# -----------------------------------------------------------------------------
# parameters of the solvers
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
                         "error_on_nonconvergence": False}

# Geometry paramaters
L, H, W = 1.0, 0.05, 0.02
N       = userpar["meshsize"]
hsize   = float(L/N)

# Loading Parameters
ut      = 1.0 # reference value for the loading (imposed displacement)

# Numerical parameters of the alternate minimization
maxiteration = 2000
AM_tolerance = 1e-4

modelname = "Taylor-Hood-FEM"
meshname  = modelname+"-mesh.xdmf"
simulation_params = "L_%.1f_H_%.1f_W_%.2f_h_%.4f" % (L, H, W, hsize)
savedir   = "output/"+modelname+"/"+simulation_params+"/"

if MPI.rank(MPI.comm_world) == 0:
    if os.path.isdir(savedir):
        shutil.rmtree(savedir)

# Mesh generation
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(L, H, W), N, int(N*H/L), int(N*W/L))
geo_mesh  = XDMFFile(MPI.comm_world, savedir+meshname)
geo_mesh.write(mesh)
"""
#read mesh
mesh = Mesh()
XDMF = XDMFFile(MPI.comm_world, "traction_3Dbar.xdmf")
XDMF.read(mesh)
"""
mesh.init()
ndim = mesh.geometry().dim()  # get number of space dimensions
if MPI.rank(MPI.comm_world) == 0:
    print ("the dimension of mesh: {0:2d}".format(ndim))

# regularization paramater
ell = Constant(5.0*hsize) # damage paramaters

# Bulk modulus
kappa = 2.0*mu*(1.0+nu)/(3.0*(1.0-2.0*nu))
#-----------------------------------------------------------------------------
p0 = -(3.0*8.0/3.0*mu*5.0*hsize+Gc)/(8.0/3.0*mu*5.0*hsize)
q0 = 2.0
tc = 2.*sqrt(-p0/3.0)*cos(1./3.*acos(3.0*q0/2.0/p0*sqrt(-3.0/p0)))-1.0

if MPI.rank(MPI.comm_world) == 0:
  print("The critical loading: [{}]".format(tc))
  print("The kappa/mu: {0:4e}".format(kappa/mu))
  print("The mu/Gc: {0:4e}".format(mu/Gc))

# ----------------------------------------------------------------------------
# Define boundary sets for boundary conditions
# ----------------------------------------------------------------------------
def left_boundary(x, on_boundary):
    return on_boundary and near(x[0], 0.0, 0.1 * hsize)

def right_boundary(x, on_boundary):
    return on_boundary and near(x[0], L, 0.1 * hsize)
"""
## when using "pointwise", the boolean argument on_boundary
## in SubDomain::inside will always be false
def left_pinponts1(x, on_boundary):
    return  (x[0]-0.)**2 + (x[1]-H)**2 + (x[2]-0.)**2 < (0.1*mesh.hmax())**2
def right_pinponts1(x, on_boundary):
    return  (x[0]-L)**2 + (x[1]-H)**2 + (x[2]-0.)**2 < (0.1*mesh.hmax())**2
def left_pinponts2(x, on_boundary):
    return  (x[0]-0.)**2 + (x[1]-0.)**2  < (0.1*mesh.hmax())**2
def right_pinponts2(x, on_boundary):
    return  (x[0]-L)**2 + (x[1]-0.)**2  < (0.1*mesh.hmax())**2

for x in mesh.coordinates():
  if (x[0]-0.)**2 + (x[1]-H)**2 + (x[2]-0.)**2 < (0.1*mesh.hmax())**2:
      print('%s is on x = L' % x)
  if (x[0]-L)**2 + (x[1]-H)**2 + (x[2]-0.)**2 < (0.1*mesh.hmax())**2:
      print('%s is on x = R' % x)
  if (x[0]-0.)**2 + (x[1]-0.)**2  < (0.1*mesh.hmax())**2:
      print('%s is on x = L3' % x)
  if (x[0]-L)**2 + (x[1]-0.)**2  < (0.1*mesh.hmax())**2:
      print('%s is on x = R3' % x)
"""
# ----------------------------------------------------------------------------
# Constitutive functions of the damage model
# ----------------------------------------------------------------------------
# Constitutive functions of the damage model
def w(alpha):
    return alpha

def a(alpha):
    return (1.0-alpha)**2

def b(alpha):
    return (1.0-alpha)**3

# ----------------------------------------------------------------------------
# Variational formulation
# ----------------------------------------------------------------------------
# Create function space for elasticity + Damage
# Taylor-Hood space for incompressible elasticity
P1      = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
P2      = VectorElement("Lagrange", mesh.ufl_cell(), 2)
TH      = MixedElement([P2,P1])
V_u     = FunctionSpace(mesh, TH)

V_alpha = FunctionSpace(mesh, "Lagrange", 1)

# Define the function, test and trial fields
w_p     = Function(V_u)
u_p     = TrialFunction(V_u)
v_q     = TestFunction(V_u)
(u, p)  = split(w_p)
alpha   = Function(V_alpha)
dalpha  = TrialFunction(V_alpha)
beta    = TestFunction(V_alpha)
"""
# --------------------------------------------------------------------
# Dirichlet boundary condition
# --------------------------------------------------------------------
u0 = Expression("0.0", degree=0)
u1 = Expression( "t",  t=0.0, degree=0)
# bc - u (imposed displacement)
bc_u0 = DirichletBC(V_u.sub(0).sub(0), u0, left_boundary)
bc_u1 = DirichletBC(V_u.sub(0).sub(0), u1, right_boundary)
bc_u2 = DirichletBC(V_u.sub(0).sub(1), u0, left_pinponts2,  method="pointwise")
bc_u3 = DirichletBC(V_u.sub(0).sub(1), u0, right_pinponts2, method="pointwise")
bc_u4 = DirichletBC(V_u.sub(0).sub(2), u0, left_pinponts1,  method="pointwise")
bc_u5 = DirichletBC(V_u.sub(0).sub(2), u0, right_pinponts1, method="pointwise")
bc_u  = [bc_u0, bc_u1, bc_u2, bc_u3, bc_u4, bc_u5]
"""
u0 = Expression(["0.0","0.0","0.0"], degree=0)
u1 = Expression("t",  t=0.0, degree=0)
# bc - u (imposed displacement)
bc_u0 = DirichletBC(V_u.sub(0), u0, left_boundary)
bc_u1 = DirichletBC(V_u.sub(0).sub(0), u1, right_boundary)
bc_u = [bc_u0, bc_u1]

# bc - alpha (zero damage)
bc_alpha0 = DirichletBC(V_alpha, 0.0, left_boundary)
bc_alpha1 = DirichletBC(V_alpha, 0.0, right_boundary)
bc_alpha = [bc_alpha0, bc_alpha1]

# --------------------------------------------------------------------
# Define the energy functional of damage problem
# --------------------------------------------------------------------
# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
C = F.T*F                   # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C) # I1
I2 = ((tr(C))**2-tr(C*C))/2.
J  = det(F) # I3 = det(C)=J**2

# Define some parameters for the eigenvalues
d_par = tr(C)/3.
e_par = sqrt(tr((C-d_par*I)*(C-d_par*I))/6.)
f_par = conditional(eq(tr((C-d_par*I)*(C-d_par*I)),0.), 0.*I, (1./e_par)*(C-d_par*I))

# Bound the argument of 'acos ' both from above and from below
# very important
g_par0 = det(f_par)/2.
g_par1 = conditional( ge(g_par0,  1.-DOLFIN_EPS),  1.-DOLFIN_EPS, g_par0 )
g_par  = conditional( le(g_par0, -1.+DOLFIN_EPS), -1.+DOLFIN_EPS, g_par1 )
h_par = acos(g_par)/3.

# Define the eigenvalues of C (principal stretches)
# lmbda1s <= lmbda2s<= lmbda3s
lmbda1s = d_par+2.*e_par*cos(h_par+2.*pi/3.)
lmbda2s = d_par+2.*e_par*cos(h_par+4.*pi/3.)
lmbda3s = d_par+2.*e_par*cos(h_par+6.*pi/3.)

# Define the Heaviside step function
# H(x) = 1, x>= 0; H(x) = 0, x< 0;
def Heaviside(x):
    return conditional(ge(x, 0.), 1., 0.)

# Elastic energy decomposition
W_act = (a(alpha)+k_ell)*Heaviside(lmbda1s-1.)*(mu/2.)*(lmbda1s-1.-ln(lmbda1s))*dx \
      + (a(alpha)+k_ell)*Heaviside(lmbda2s-1.)*(mu/2.)*(lmbda2s-1.-ln(lmbda2s))*dx \
      + (a(alpha)+k_ell)*Heaviside(lmbda3s-1.)*(mu/2.)*(lmbda3s-1.-ln(lmbda3s))*dx \
      - b(alpha)*p*ln(J)*Heaviside(J-1.)*dx - 1./(2.*kappa)*p**2*Heaviside(J-1.)*dx

W_pas = (mu/2.)*Heaviside(1.-lmbda1s)*(lmbda1s-1.-ln(lmbda1s))*dx \
      + (mu/2.)*Heaviside(1.-lmbda2s)*(lmbda2s-1.-ln(lmbda2s))*dx \
      + (mu/2.)*Heaviside(1.-lmbda3s)*(lmbda3s-1.-ln(lmbda3s))*dx \
      - p*ln(J)*Heaviside(1.-J)*dx - 1./(2.*kappa)*p**2*Heaviside(1.-J)*dx

elastic_energy = (a(alpha)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx-b(alpha)*p*ln(J)*dx-1./(2.*kappa)*p**2*dx

elastic_energy2= (a(alpha)+k_ell)*(mu/2.0)*(lmbda1s+lmbda2s+lmbda3s-3.0-ln(lmbda1s*lmbda2s*lmbda3s))*dx-b(alpha)*p*ln(J)*dx-1./(2.*kappa)*p**2*dx

body_force        = Constant((0., 0., 0.))
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

# --------------------------------------------------------------------
# Define the energy functional of damage problem
# --------------------------------------------------------------------
alpha_0 = interpolate(Expression("0.", degree=0), V_alpha)  # initial (known) alpha
z = sympy.Symbol("z", positive=True)
c_w = float(4 * sympy.integrate(sympy.sqrt(w(z)), (z, 0, 1)))
dissipated_energy = Gc/float(c_w)*(w(alpha)/ell + ell*dot(grad(alpha), grad(alpha)))*dx
damage_functional = W_act + W_pas + dissipated_energy

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

# Set up the solvers
solver_alpha  = PETScTAOSolver()
solver_alpha.parameters.update(tao_solver_parameters)
# info(solver_alpha.parameters,True) # uncomment to see available parameters

#alpha_lb = interpolate(Expression("x[0]<=0.5*L & near(x[1], 0.0, tol) ? 1.0 : 0.0", \
#                                  degree=0, L= L, tol=1.2*cra_w), V_alpha)  # initial (known) alpha
alpha_lb = interpolate(Expression("0.", degree=0), V_alpha)  # lower bound, set to 0
alpha_ub = interpolate(Expression("1.", degree=0), V_alpha)  # upper bound, set to 1

# loading and initialization of vectors to store time datas
load_multipliers  = np.linspace(userpar["load_min"], userpar["load_max"], userpar["load_steps"])
energies          = np.zeros((len(load_multipliers), 5))
iterations        = np.zeros((len(load_multipliers), 2))

# set the saved data file name
(u, p)      = w_p.split()
file_u      = XDMFFile(MPI.comm_world, savedir + "/u.xdmf")
file_u.parameters["rewrite_function_mesh"]     = False
file_u.parameters["flush_output"]              = True
file_p      = XDMFFile(MPI.comm_world, savedir + "/p.xdmf")
file_p.parameters["rewrite_function_mesh"]     = False
file_p.parameters["flush_output"]              = True
file_alpha  = XDMFFile(MPI.comm_world, savedir + "/alpha.xdmf")
file_alpha.parameters["rewrite_function_mesh"] = False
file_alpha.parameters["flush_output"] = True
# write the parameters to file
File(savedir+"/parameters.xml") << userpar

# ----------------------------------------------------------------------------
# Solving at each timestep
# ----------------------------------------------------------------------------
for (i_t, t) in enumerate(load_multipliers):
    u1.t = t * ut
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
        #u_file = File("displacement.pvd")
        #u_file << u
        #a_file = File("damage.pvd")
        #a_file << alpha
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
    alpha.rename("Damage", "alpha")
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")

    # Dump solution to file
    file_alpha.write(alpha, t)
    file_u.write(u, t)
    file_p.write(p, t)

    # ----------------------------------------
    # Some post-processing
    # ----------------------------------------
    # Save number of iterations for the time step
    iterations[i_t] = np.array([t, iteration])

    # Calculate the energies
    elastic_energy_value = assemble(elastic_energy)
    elastic_energy_value2 = assemble(elastic_energy2)
    surface_energy_value = assemble(dissipated_energy)
    volume_ratio         = assemble(J/(L*H*W)*dx)

    energies[i_t] = np.array([t, elastic_energy_value, surface_energy_value, elastic_energy_value+\
    	                      surface_energy_value, volume_ratio])

    if MPI.rank(MPI.comm_world) == 0:
        print("\nEnd of timestep {0:3d} with load multiplier {1:4f}".format(i_t, t))
        print("\nElastic and Surface Energies: [{0:6f},{1:6f}]".format(elastic_energy_value, surface_energy_value))
        print("\nElastic and Surface Energies: [{},{}]".format(elastic_energy_value, surface_energy_value))
        print("-----------------------------------------")
        print("\nElastic Energies [1,2]: [{},{}]".format(elastic_energy_value, elastic_energy_value2))
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
