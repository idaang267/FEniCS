# Small strain, linear elastic solution for a 2D isotropic linear elastic medium
# in plane stress or plane strain, in a traditional displacment-based
# finite element formulation

from dolfin import *
import numpy as np

# Parameters
# -----------------------------------------------------------------------------
# Set the user parameters
parameters.parse()
userpar = Parameters("user")
userpar.add("model", "plane_stress")
userpar.add("MaxStep", 50)
userpar.add("MaxRho", 0.2)
# Parse command-line options
userpar.parse()

# Chose either plane strain or plane stress
model = userpar["model"]

# Cantiliver beam modeled as a 2D medium of dimensions (L x H)
L, H = 10., 1.
# Mesh density parameters
Nx, Ny = 100, 10
# Body force
MaxRho = userpar["MaxRho"]
MaxStep = userpar["MaxStep"]
ArrRho = np.linspace(0, MaxRho, MaxStep)
f = Expression(("0", "-Rho"), Rho=0.0, degree=1)

# Material Parameters
E = Constant(1000)               # Young's Modulus
nu = Constant(0.3)              # Poisson's Ratio
mu = E/(2*(1+nu))               # Lamé Coefficient
lmbda = E*nu/((1+nu)*(1-2*nu))  # Lamé Coefficient

# Consider a 2D model in either plane strain or plane stress conditions. Either
# condition means we work with a 2D displacement vecotr u = (u_x, u_y)
if model == "plane_strain":
    lmbda = lmbda
else:
    lmbda = 2*mu*lmbda/(lmbda+2*mu)

SimPar = "S_%.0d" % (MaxStep)
ResDir = "../Result/" + model + '/' + SimPar + '/'

# Mesh density defined with a criss-crossed structured mesh
mesh = RectangleMesh(Point(0, 0), Point(L, H), Nx, Ny, "crossed")

# Variational Problem
# ---------------------
# Define the function space
V_CG2 = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
T_DG0 = TensorFunctionSpace(mesh, "DG", degree=0)
v = TestFunction(V_CG2)                     # Test Function
u = Function(V_CG2, name="Displacement")    # Function named for paraview output
sig = Function(T_DG0, name="Stress")

# Boundary Condition
# ---------------------
def left(x, on_boundary):
    return near(x[0],0.)

bc = DirichletBC(V_CG2, Constant((0.,0.)), left)

# Constitutive relationships
# ---------------------
# Strain operator computes the 2x2 plane components of the symmetrized gradient
# tensor of any 2D vectorial field
def eps(v):
    return sym(grad(v))
# Stress operator computed depending on plane strain or stress
def sigma(v):
    return lmbda*tr(eps(v))*Identity(2) + 2.0*mu*eps(v)

# Solve for displacement
F = inner(sigma(u),eps(v))*dx - inner(f,v)*dx

# Fields can be exported in a suitable format for visualization in Paraview
file_results = XDMFFile(ResDir + "/Results.xdmf")
# Control certain parameters
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

# Postprocessing
PostTxt = np.zeros((MaxStep, 4))

for (Step, Rho) in enumerate(ArrRho):
    # Print outs to track code progress
    print("\033[1;32m--- Step {0:2d}: Rho = {1:2f} ---\033[1;m".format(Step, Rho))

    solve(F == 0, u, bc)

    # Projecting the computed stress onto the solution from beam theory in order to
    # evaluate pointwise values
    sig.assign(project(sigma(u), T_DG0))

    # Only one time step, 0
    file_results.write(u, Step)
    file_results.write(sig, Step)

    # Postprocessing:
    PostTxt[Step, :] = np.array([Step, Rho, -u(L,H/2.)[1], 3*Rho*L**4/2/E/H**3])
    np.savetxt(ResDir + '/PostProc.txt', PostTxt)

    f.Rho = Rho     # Update expression class

# # Validate the results
# print("Maximal deflection:", -u(L,H/2.)[1])
# print("Beam theory deflection:", float(3*MaxRho*L**4/2/E/H**3))
# print("Stress at (0,H):", sig(0, H))
