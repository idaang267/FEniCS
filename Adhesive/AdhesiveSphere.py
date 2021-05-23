# Hyperelasticity

# Import modules
from dolfin import *
import matplotlib.pyplot as plt     # For visualization
import numpy as np

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"
# The following two fields were in the original hyperelasticity demo
parameters["form_compiler"]["cpp_optimize"] = False
parameters["form_compiler"]["representation"] = "uflacs"
# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          'absolute_tolerance':1e-5,
                                          'relative_tolerance':1e-5,
                                          "maximum_iterations": 20,
                                          "report": False,
                                          "error_on_nonconvergence": True}}

solver_par = NonlinearVariationalSolver.default_parameters()
solver_par.rename("solver")

# Control the number of loop iterations through user parameters
parameters.parse()
user_par = Parameters("user")
user_par.add("u_min",0.)        # Displacement minimum
user_par.add("u_max",0.5)       # Displacement maximum
user_par.add("u_nsteps",5)     # Displacement step number
# Parse command-line options
user_par.parse()

# Elasticity parameters
E, nu = 10.0, 0.3
mu = Constant(E/(2*(1 + nu)))
lmbda = Constant(E*nu/((1+nu)*(1-2*nu)))

r_d = 7.1
t = 0.2
hsize = 0.4

# Naming parameters for saving output
modelname = "HydrogelPlug"
meshname  = modelname + "-mesh.xdmf"
simulation_params = "Test"
savedir   = "output/"+"/"+simulation_params+"/"

# Define mesh and mixed function space
# mesh = UnitCubeMesh(20, 10, 10)
mesh = Mesh("SphereWPlug.xml")
physical = MeshFunction("size_t", mesh, "SphereWPlug_physical_region.xml")
facet = MeshFunction("size_t", mesh, "SphereWPlug_facet_region.xml")

class Inner(SubDomain):
    def inside(self, x, on_boundary):
        r_sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
        return r_sq <= r_d*r_d and on_boundary

class pin_point_b(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0, 0.4*hsize) and near(x[1], 0, 0.4*hsize) and near(x[2], -r_d, 0.4*hsize)

class pin_point_s(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], r_d, 0.4*hsize) and near(x[1], 0, 0.4*hsize) and near(x[2], 0, 0.4*hsize)

Inner = Inner()
pin_point_b = pin_point_b()
pin_point_s = pin_point_s()


lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 3)

# show lines of interest
lines.set_all(0)
Inner.mark(lines, 1)
file_results = XDMFFile(savedir + "/" + "lines.xdmf")
file_results.write(lines)

points.set_all(0)
pin_point_b.mark(points, 1)
pin_point_s.mark(points, 1)
file_results = XDMFFile(savedir + "/" + "points.xdmf")
file_results.write(points)

# show boundary  of interest
geo_mesh = XDMFFile(MPI.comm_world, savedir + meshname)
geo_mesh.write(mesh)

geo_mesh = XDMFFile(MPI.comm_world, savedir + "physical.xdmf")
geo_mesh.write(physical)

geo_mesh = XDMFFile(MPI.comm_world, savedir + "facet.xdmf")
geo_mesh.write(facet)

# Taylor-Hood space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)  # Displacement
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)  # Hydrostatic Pressure
TH = MixedElement([P2,P1])
V  = FunctionSpace(mesh, TH)

# Define trial and test functions and
# unknown solutions (u,p) = (displacement, hydrostatic pressure)
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
w  = Function(V)                 # Solutions (u,p) from previous iteration

# NOTE that this is a shallow copy not a deep copy
(v_u, v_p) = split(v)            # Split the test functions
(u, p) = split(w)                # Split the solutions to (u,p)

# Define Dirichlet boundary (x = 0 or x = 1)
u00 = Expression("0.0", degree=0)
u0 = Expression(["t*x[0]/r_d","t*x[1]/r_d","t*x[2]/r_d"], t=0.0,r_d = r_d, degree=0)
c = Constant((0.0, 0.0, 0.0))
# The Dirichlet BCs are specified in a subspace
bc_in1 = DirichletBC(V.sub(0), u0, facet, 1)
bc_in2 = DirichletBC(V.sub(0), u0, facet, 2)
bc_px1 = DirichletBC(V.sub(0).sub(0), u00, pin_point_b)
bc_py1 = DirichletBC(V.sub(0).sub(1), u00, pin_point_b)
bc_py2 = DirichletBC(V.sub(0).sub(1), u00, pin_point_s)
bc_pz2 = DirichletBC(V.sub(0).sub(2), u00, pin_point_s)
bcs = [bc_in1, bc_in2, bc_px1, bc_py1, bc_py2, bc_pz2]

# Definitions of constants from hyperelasticity demo
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary

# Kinematics
d = len(u)                      # Spatial dimension
I = Identity(d)                 # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Where P = dW/dF:
def P(u):
    return mu*F + p*J*inv(F.T)

# Specify the quadrature degree for efficiency
F = (inner(P(u), grad(v_u)) + inner(J-1., v_p))*dx(metadata={"quadrature_degree": 4}) \
    - dot(B, v_u)*dx - dot(T, v_u)*ds

# Directional derivative
J_o = derivative(F, w, du)

# Solve variational problem
varproblem = NonlinearVariationalProblem(F, w, bcs, J=J_o)
solver = NonlinearVariationalSolver(varproblem)
solver.parameters.update(snes_solver_parameters)
info(solver.parameters, False)
(iter, converged) = solver.solve()

# Loading parameter (list of values of surface tension for the simulations)
u_list = np.linspace(user_par["u_min"],user_par["u_max"],user_par["u_nsteps"])

# Save results to an .xdmf file since we have multiple fields
file_results = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")

# Solve with Newton solver for each displacement value using the previous
# solution as a starting point
for t in u_list:
    # Update the displacement value
    u0.t = t
    # Solve the nonlinear problem (using Newton solver)
    solver.solve()

    # Note that this is now a deep copy not a shallow copy
    (u, p) = w.split()
    u.rename("Displacement", "u")
    p.rename("Pressure","p")
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    file_results.write(u,t)
    file_results.write(p,t)
