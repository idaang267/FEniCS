# Hertzian Contact Problem with Penalty Formulation
# -----------------------------------------------------------------------------
from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt

# # Iterative Conjugate-Gradient (cg) solver is chosen for the Newton solver
# parameters["newton_solver"]["linear_solver"] = "cg"
# Incomplete Lower-Upper (ilu) solver chosen for the Newton solver
# parameters["newton_solver"]["preconditioner"] = "ilu"

# PETSc SNES solver: non-linear solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          'absolute_tolerance':1e-5,
                                          'relative_tolerance':1e-5,
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"
# The following two fields were in the original hyperelasticity demo
parameters["form_compiler"]["cpp_optimize"] = True
# For processing UFLACS - unified form language expressions
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["allow_extrapolation"] = True

solver_par = NonlinearVariationalSolver.default_parameters()
solver_par.rename("solver")

# Element-wise projection using LocalSolver
def local_project(v, V, u=None):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

# This is the top domain that the indenter will contact
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0.5) and on_boundary

# Bottom which will be fully fixed
class Bot(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], -0.5) and on_boundary

# Define symmetry conditions on the x = 0 and y = 0 surfaces
class SymX(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], -0.5) and on_boundary

class SymY(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.5) and on_boundary

Top, Bot = Top(), Bot()
SymX, SymY = SymX(), SymY()

# Parameters
# -----------------------------------------------------------------------------
N, N_plane = 10, 20                              # Mesh refinement
TotSteps = 15                   # Total number of time steps for indentation
# Material Parameters
E = Constant(10.)                   # Young's Modulus
nu = Constant(0.3)                  # Poisson's Ratio
mu = E/(2*(1+nu))                   # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))      # Lamé Coefficient
# A large penalty parameter deteriorates the problem conditioning so the
# solving time will drastically increase and the problem can fail
pen = Constant(1e4)
# Define the indenter (obstacle) with a parabola
R = 0.25
DepthExp = np.linspace(0, 0.2, TotSteps)
# Obstacle in the X and Y plane
obstacle = Expression("-d+(pow(x[0],2)+pow(x[1],2))/(2*R)", d=0., R=R, degree=2)
# Postprocessing
PointInt = 4*N_plane + 1  # Number of points across the profile to examine (odd not even)
PlotMod = 2         # Controls the profiles to plot

# File output name
ResDir = "../Result/3DPenalty/"
name = ResDir + "Penalty.xdmf"

# Define mesh and function spaces
# -----------------------------------------------------------------------------
# z is in the vertical direction, the middle of the contact area is x=y=z = 0

# N//2 is python syntax for floor division in case N is not even
# mesh = UnitCubeMesh.create(N, N, N//2, CellType.Type.hexahedron)
mesh = BoxMesh(Point(-.5, -.5, -.5), Point(.5, .5, .5), N_plane, N_plane, N)
# Exterior facets
subdomains = MeshFunction("size_t", mesh, 3, mesh.domains())
facets = MeshFunction("size_t", mesh, 2)

# All facets are initially marked as subdomain 0 on all surfaces of the cube
facets.set_all(0)
# Top exterior facets are marked as 1 using the top class defined above
Top.mark(facets, 1)
Bot.mark(facets, 2)
SymX.mark(facets, 3)
SymY.mark(facets, 4)
XDMFFile(ResDir + "boundaries.xdmf").write(facets)

# Measure redefines ds
ds = Measure('ds', subdomain_data=facets)

# Variational Problem
# -----------------------------------------------------------------------------
# Define the function spaces for displacement, gap, and contact pressure
V_CG1 = VectorFunctionSpace(mesh, "CG", 1)      # Displacement
CG1 = FunctionSpace(mesh, "CG", 1)           # Gap
# Discontinuous Galerkin (DG) for discontinuous basis functions for mesh and
# space adaptivity
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)
DG0 = FunctionSpace(mesh, "DG", 0)           # Contact Pressure

# Define the trial and test functions
du = TrialFunction(V_CG1)
v = TestFunction(V_CG1)
# Define the three functions we will output (names for Paraview)
u = Function(V_CG1, name="Displacement")
gap = Function(CG1, name="Gap")
p = Function(DG0, name="Contact pressure")
CurStress = Function(T_DG0, name="Stress")

# Boundary Conditions
# -----------------------------------------------------------------------------
bc_b = DirichletBC(V_CG1, Constant((0.0, 0.0, 0.0)), Bot)
# Subclass the displacement domain to set x = 0
bc_x = DirichletBC(V_CG1.sub(0), Constant(0.0), SymX)
# Subclass the displacement domain to set y = 0
bc_y = DirichletBC(V_CG1.sub(1), Constant(0.0), SymY)
bc = [bc_b, bc_x, bc_y]

# Constitutive relation
# -----------------------------------------------------------------------------
def eps(v):     # Strain function
    return sym(grad(v))
def sigma(v):   # Stress function defined using the Neo-Hookean Model
    return lmbda*tr(eps(v))*Identity(3) + 2.0*mu*eps(v)
def ppos(x):    # Definition of The Mackauley bracket <x>_+
    return (x+abs(x))/2.
def gap_f(u):     # Gap function
    return u[2] - obstacle

# Weak form and solver type
# -----------------------------------------------------------------------------
# Weak form with penalty term added multiplied by parameter
WF = inner(sigma(u), eps(v))*dx + pen*dot(v[2], ppos(gap_f(u)))*ds(1)
# Jacobian
Ju = derivative(WF, u, du)

# Define a non-linear solver
problem = NonlinearVariationalProblem(WF, u, bc, J=Ju)
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(snes_solver_parameters)

# Fields can be exported in a suitable format for visualization in Paraview
file = XDMFFile(name)
# Control parameters: If file crashes, this will save output
file.parameters["flush_output"] = True
# Makes sure each function (u, gap, p) we are saving share the same mesh
file.parameters["functions_share_mesh"] = True

# Arrays for storing data
Points = np.linspace(-0.5, 0.5, PointInt)       # Points along the profile
TempArray = np.zeros((1, len(Points)))          # Array for storing displacements along profile
Profile = np.zeros((TotSteps+1, len(Points)))  # Saving profile
Profile[0,:] = Points[:]

# Solution
# -----------------------------------------------------------------------------
for (step, depth) in enumerate(DepthExp):
    # Print outs to track code progress
    print("\033[1;32m--- Step {0:2d}: depth = {1:2f} ---\033[1;m".format(step, depth))

    solver.solve()      # solve for each step

    # Compute the pressure and gap by projecting onto the respective function space
    # local_project(-sigma(u)[2, 2], DG0, p)
    local_project(pen*u[2], DG0, p)              # Pressure
    local_project(ppos(gap_f(u)), CG1, gap)      # Gap
    local_project(sigma(u), T_DG0, CurStress)    # Current stress

    # Only one time step, 0, but this still counts as multiple fields
    file.write(u, step)
    file.write(gap, step)
    file.write(p, step)

    # Postprocessing: int is the interval along the profile
    for (i, int) in enumerate(Points):
        TempArray[:, i] = u(int,0.5,0)[2]
    Profile[step+1,:] = TempArray[:]
    # Save for future post-processing
    np.savetxt(ResDir + '/Profile.txt', Profile)

    # Update depth in expression class
    obstacle.d = depth

# Plot profile evolution
plt.figure()
ii = 1
while ii <= TotSteps:
    if ii%PlotMod == 0:
        plt.plot(Points, Profile[ii,:], label='Step {0:2d}'.format(ii))
    ii += 1

plt.legend(loc='best', frameon=False)
plt.xlabel("Steps")
plt.ylabel("Displacement")
plt.savefig(ResDir + "/ProfileEvolution.pdf", transparent=True)
plt.close()

# # Analytical solution for comparison
# a = sqrt(R*d)
# F = 4/3.*float(E)/(1-float(nu)**2)*a*d
# p0 = 3*F/(2*pi*a**2)
# print("Maximum pressure FE: {0:8.5f} Hertz: {1:8.5f}".format(max(np.abs(p.vector().get_local())), p0))
# print("Applied force FE: {0:8.5f} Hertz: {1:8.5f}".format(4*assemble(p*ds(1)), F))
