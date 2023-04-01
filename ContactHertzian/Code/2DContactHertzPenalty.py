# Hertzian Contact Problem with Penalty Formulation
# -----------------------------------------------------------------------------
from dolfin import *
# from multiphenics import *
from mshr import *          # For generating domains
import numpy as np
import matplotlib.pyplot as plt

parameters["ghost_mode"] = "shared_facet" # required by dS

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
        return near(x[1], 1.0) and on_boundary

# bottom = CompiledSubDomain("near(x[1], 0.0) && on_boundary")
class Bot(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

Top, Bot = Top(), Bot()

# Parameters
# -----------------------------------------------------------------------------
# Set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("N", 30)
userpar.add("TotSteps", 10)
userpar.add("pen", 1e4)
# Parse command-line options
userpar.parse()

N = userpar["N"]                # Mesh Resolution
TotSteps = userpar["TotSteps"]  # Total number of time steps for indentation
# A large penalty parameter deteriorates the problem conditioning so the
# solving time will drastically increase and the problem can fail
pen = userpar["pen"]
# Material Parameters
E = Constant(10.0)                  # Young's Modulus
nu = Constant(0.3)                  # Poisson's Ratio
mu = E/(2*(1+nu))                   # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))      # Lamé Coefficient
# Define the indenter (obstacle) with a parabola
R = 0.25                                    # Radius
DepthExp = np.linspace(0, 0.2, TotSteps)    # Max depth can be specified here
# Indenter
obstacle = Expression("-depth+(pow(x[0]-0.5,2))/(4*R)", depth=0.0, R=R, degree=2)
# Postprocessing
PointInt = 4*N + 1  # Number of points across the profile to examine (odd not even)
PlotMod = 2         # Controls the profiles to plot

# File output name
SimPar = "Pen_%.0e" % (pen)
ResDir = "../Result/2DPenalty/" + SimPar + '/'
name = ResDir + "Penalty.xdmf"

# Define mesh and function spaces
# -----------------------------------------------------------------------------
# z is in the vertical direction, the middle of the contact area is x=y=z = 0
domain = Rectangle(Point(0.0, 0.0), Point(1.0, 1.0))
mesh = generate_mesh(domain, N)
# Create subdomains and boundaries
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

# Initially mark boundaries as 0 on all surfaces of the cube
boundaries.set_all(0)
# Mark boundaries using the classes defined above
Top.mark(boundaries, 1)     # mark top
Bot.mark(boundaries, 2)     # Mark bottom
XDMFFile(ResDir + "boundaries.xdmf").write(boundaries)

# Measure redefines ds
dx = dx(subdomain_data=subdomains)
ds = ds(subdomain_data=boundaries)

# Define the function spaces for displacement, gap, and contact pressure
V_CG1 = VectorFunctionSpace(mesh, "CG", 1)      # Displacement
# Discontinuous Galerkin (DG) for discontinuous basis functions for mesh and space adaptivity
DG0 = FunctionSpace(mesh, "DG", 0)           # Contact Pressure
T_DG0 = TensorFunctionSpace(mesh, "DG", 0)

# Variational Problem
# -----------------------------------------------------------------------------
# Define the trial and test functions
du = TrialFunction(V_CG1)
v = TestFunction(V_CG1)
u = Function(V_CG1, name="Displacement")        # Displacement
p = Function(DG0, name="Pressure")              # Pressure
gap = Function(DG0, name="gap")                 # Gap
CurStress = Function(T_DG0, name="Stress")      # Cauchy (Current) Stress

# Boundary Conditions
# -----------------------------------------------------------------------------
# Fix bottom
bc = DirichletBC(V_CG1, Constant((0.0, 0.0)), Bot)

# Constitutive relation
# -----------------------------------------------------------------------------
def eps(v):     # Strain function
    return sym(grad(v))
def sigma(v):   # Stress function defined using the Neo-Hookean Model
    return lmbda*tr(eps(v))*Identity(2) + 2.0*mu*eps(v)
def ppos(x):    # Definition of The Mackauley bracket <x>_+
    return (x+abs(x))/2.
def gap_f(u):
    return u[1]-obstacle

# Weak form and solver type
# -----------------------------------------------------------------------------
# Weak form with penalty term added multiplied by parameter
WF = inner(sigma(u), eps(v))*dx + Constant(pen)*dot(v[1],ppos(gap_f(u)))*ds(1)
# Jacobian
Ju = derivative(WF, u, du)

# Define a non-linear solver
problem = NonlinearVariationalProblem(WF, u, bc, J=Ju)
solver = NonlinearVariationalSolver(problem)
# Iterative Conjugate-Gradient (cg) solver is chosen for the Newton solver
solver.parameters["newton_solver"]["linear_solver"] = "cg"
# Incomplete Lower-Upper (ilu) solver chosen for the Newton solver
solver.parameters["newton_solver"]["preconditioner"] = "ilu"

# Solution
# -----------------------------------------------------------------------------
# Fields can be exported in a suitable format for visualization in Paraview
file = XDMFFile(name)                           # .xdmf file for multiple fields
file.parameters["flush_output"] = True
file.parameters["functions_share_mesh"] = True

# Arrays for storing data
Points = np.linspace(0, 1.0, PointInt)         # Points along the profile
DispArray = np.zeros((1, len(Points)))         # Array for storing displacements along profile
GapArray = np.zeros((1, len(Points)))
DispProf = np.zeros((TotSteps+1, len(Points)))  # Saving profile
GapProf = np.zeros((TotSteps+1, len(Points)))
PostTxt = np.zeros((TotSteps+1, 2))            # Save displacement and gap across profile

# Loop to increase indenter depth slowly
for (step, depth) in enumerate(DepthExp):
    # Print outs to track code progress
    print("\033[1;32m--- Step {0:2d}: depth = {1:2f} ---\033[1;m".format(step, depth))

    solver.solve()              # Solve for each step

    # Project pressure
    local_project(pen*u[1], DG0, p)
    local_project(ppos(gap_f(u)), DG0, gap)
    local_project(sigma(u), T_DG0, CurStress)

    # Write to file
    file.write(u, step)
    file.write(p, step)
    file.write(gap, step)
    file.write(CurStress, step)

    # Postprocessing: int is the interval along the profile
    PostTxt[step, :] = np.array([step, depth])

    for (i, int) in enumerate(Points):
        DispArray[:, i] = u(int, 1)[1]
        GapArray[:, i] = gap(int,1)

    DispProf[step+1,:] = DispArray[:]
    GapProf[step+1,:] = GapArray[:]

    # Save for future post-processing
    np.savetxt(ResDir + '/DispProf.txt', DispProf)
    np.savetxt(ResDir + '/GapProf.txt', GapProf)
    np.savetxt(ResDir + '/PostProc.txt', PostTxt)

    # Update depth in expression class
    obstacle.depth = depth

# # Plot profile evolution
# plt.figure()
# ii = 1
# while ii <= TotSteps:
#     if ii%PlotMod == 0:
#         plt.plot(Points, Profile[ii,:], label='Step {0:2d}'.format(ii))
#     ii += 1
#
# plt.legend(loc='best', frameon=False)
# plt.xlabel("Steps")
# plt.ylabel("Displacement")
# plt.savefig(ResDir + "/ProfileEvolution.pdf", transparent=True)
# plt.close()
