# This demo program solves a hyperelastic problem with surface energy to model
# the large deformations of a capillary solids
#-------------------------------------------------------------
# Author: Corrado Maurini (corrado.maurini@upmc.fr)
# First version: 08/2011
# last modified 08/2013

# Load the required modules
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt     # For visualization
import numpy as np
import os

set_log_level(INFO)

# Optimization options for the form compiler
#parameters["num_threads"] = 1
parameters["mesh_partitioner"] = "SCOTCH"
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["log_level"] = INFO
parameters["allow_extrapolation"] = True
ffc_options = {"optimize": True, \
                "eliminate_zeros": True, \
                "precompute_basis_const": True, \
                "precompute_ip_const": True, \
                "quadrature_degree": 2}

# Some user parameters
user_par = Parameters("user")
user_par.add("bounds_xmin",-0.5)
user_par.add("bounds_xmax",0.5)
user_par.add("bounds_ymin",-0.5)
user_par.add("bounds_ymax",0.5)
user_par.add("bounds_zmin",0.)
user_par.add("bounds_zmax", 2.5)
user_par.add("fe_order_u",1)
user_par.add("fe_order_p",1)

# Control the number of loop iterations
user_par.add("gamma_min",0.)
user_par.add("gamma_max",.2)
user_par.add("gamma_nsteps",10)
user_par.add("mesh_ref",10)

# Control specifications for saving files
user_par.add("save_dir","results")
user_par.add("output_type","pvd")
user_par.add("plot",True)

# Non-linear solver parameters
solver_par = NonlinearVariationalSolver.default_parameters()
solver_par.rename("solver")
solver_par["symmetric"]=True
# solver_par["linear_solver"]="umfpack"           # use "mumps" in parallel
# solver_par["lu_solver"]["same_nonzero_pattern"] = True
# solver_par["lu_solver"]["verbose"] = True
solver_par["newton_solver"]["maximum_iterations"] = 20
solver_par["newton_solver"]["relaxation_parameter"] = .8
solver_par["newton_solver"]["relative_tolerance"] = 1e-5
solver_par["newton_solver"]["absolute_tolerance"] = 1e-5

# Add user parameters in the global parameter set
parameters.add(user_par)
parameters.add(solver_par)

# Parse parameters from command line
parameters.parse()
info(parameters,True)
user_par = parameters.user

# Geometry
#-------------------------------------------------------------
# Create the geometry and the mesh
xmin,xmax = user_par.bounds_xmin,user_par.bounds_xmax
ymin,ymax = user_par.bounds_ymin,user_par.bounds_ymax
zmin,zmax = user_par.bounds_zmin,user_par.bounds_zmax
mesh_ref = user_par.mesh_ref
mesh = BoxMesh(dolfin.Point(xmin,ymin,zmin),dolfin.Point(xmax,ymax,zmax),mesh_ref,mesh_ref,mesh_ref)

# Old
#geom = Box(xmin,ymin,zmin,xmax,ymax,zmax)
#mesh = Mesh(geom,mesh_ref)

# Definition of function spaces
#-------------------------------------------------------------
# Create function space
P2 = VectorElement("CG",mesh.ufl_cell(),user_par.fe_order_u)   # Displacement
P1 = FiniteElement("CG",mesh.ufl_cell(),user_par.fe_order_p)   # Pressure
# Equal Order elements
EO = MixedElement([P2,P1])
V = FunctionSpace(mesh,EO)

# Define the sub-spaces for displacement and pressure using a hard copy
V_u = V.sub(1)
V_p = V.sub(0)

# Create test and trial functions for the variational formulation
dup = TrialFunction(V)
vq = TestFunction(V)
# Create functions to define the energy and store the results
up = Function(V)

#(p,u) = up.split()  # deep copy
(p,u) = split(up)  # shallow copy
(q,v) = TestFunctions(V)    # why is this not a shallow copy using split

# Boundary conditions
#-------------------------------------------------------------
# Mark boundary subdomains
xtol = mesh.hmin()/4.
class ClampedBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] - zmin < xtol and on_boundary
# Define boundaries with surface tension
class SurfaceBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Define the boundary conditions
zero_vector = Constant(0.0,0.0,0.0) # Used to be 3D (0.0 0.0 0.0)
clamped_boundary = ClampedBoundary()
bc1 = DirichletBC(V_u, zero_vector, clamped_boundary)
bc_u = [bc1]

# Mark facets where apply surface tension with 1
boundary_parts = MeshFunction("size_t", mesh, 0)
surface_boundary = SurfaceBoundary()
surface_boundary.mark(boundary_parts, 1)

# Redefine element of area to include informations about surface tension
ds = Measure("ds", domain = mesh, subdomain_data=boundary_parts)
#ds = ds[boundary_parts]

# Kinematics
#-------------------------------------------------------------
ndim = len(u)
#ndim = P2.cell().d      # Spatial Dimension
I = Identity(ndim)      # Identity tensor
F = I + grad(u)         # Deformation gradient
C = transpose(F)*F      # Right Cauchy-Green tensor
E = 0.5*(C - I)         # Green-Lagrange tensor

# Invariants of deformation tensors
Ic = tr(C)
J = det(F)

# Normal and tangent vectors in the reference configuration
N = FacetNormal(mesh)
# Element of area transformation operator
NansonOp = transpose(cofac(F))
# Surface element vector in the deformed configuration
deformed_N = dot(NansonOp,N)
# Norm of the surface element vector in the current configuration
current_element_of_area = sqrt(dot(deformed_N,deformed_N))

# Energy and variational formulation
#-------------------------------------------------------------
mu, lmbda = Constant(1.), Constant(1000.)   # Lamés parameters

# Bulk energy (strain energy for an almost incompressible neo-Hookean model)
bulk_energy_density = mu*(Ic - ndim) -(mu + p)*ln(J) - 1/(2*lmbda)*p**2
bulk_energy = bulk_energy_density*dx

# Surface energy
gamma = Expression("t",t=0.00)
surface_energy_density = gamma*current_element_of_area
surface_energy = surface_energy_density*ds(1)

# Total potential energy
potential_energy = bulk_energy + surface_energy

# First directional derivative of the potential energy (a linear form in the test function vq)
F = derivative(potential_energy,up,vq)

# First directional derivative of the potential energy (a bilinear form in the test function vq and the trial function dup)
dF = derivative(F,up,dup)

# Setup the variational problem
varproblem = NonlinearVariationalProblem(F, up, bc_u, J=dF,form_compiler_parameters=ffc_options)

# Set up the solver (Newton solver)
#-------------------------------------------------------------
solver = NonlinearVariationalSolver(varproblem)
solver.parameters.update(parameters.solver)

# Solve the problem
#-------------------------------------------------------------
# Loading parameter (list of values of surface tension for the simulations)
gamma_list = np.linspace(user_par.gamma_min,user_par.gamma_max,user_par.gamma_nsteps)
# Directory and files to save the results
save_dir = parameters.user.save_dir
file_u = File(save_dir+"/displacement."+parameters.user.output_type)
file_p = File(save_dir+"/pressure."+parameters.user.output_type)

# Solve with Newton solver for each value of the surface tension, using the previous solution as a starting point.
for t in gamma_list:
    # Update the value of the surface tension
    gamma.t = t
    # Solve the nonlinear problem (using Newton solver)
    solver.solve()
    # Save solution to file (readable by Paraview)
    (p,u) = up.split()
    file_u << (u,t)
    file_p << (p,t)
    # Plot and save png image
    if parameters.user.plot:
        plot_u = plot(u, mode = "displacement",title="Displacement field gamma=%.4f"%t,elevate=25.0)
        plot_u.write_png(save_dir+"/displacement_%.4f"%t)

# Save the parameters to file
File(save_dir+"/parameters.xml") << parameters

# Get timings and save to file
if MPI.process_number() == 0:
    timings_str = timings().str("Timings")
    text_file = open(save_dir+"/timings.txt", "w")
    text_file.write(timings_str)
    text_file.close()
