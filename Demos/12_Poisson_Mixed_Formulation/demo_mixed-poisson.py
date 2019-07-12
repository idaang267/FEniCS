# Mixed formulation for Poisson equation
# ======================================
#
# This demo illustrates how to solve Poisson equation using a mixed (two-field)
# formulation. In particular, it illustrates how to
#
# * Use mixed and non-continuous finite element spaces
# * Set essential boundary conditions for subspaces and H(div) spaces
# * Define a (vector-valued) expression using additional geometry information
#
# An alternative formulation of Poisson equation can be formulated by
# introducing an additional vector variable, namely the negative flux.
#
# Use the same definitions of functions and boundaries as in the demo for
# Poisson's equation.

# First, the required modules are imported
from dolfin import *
import matplotlib.pyplot as plt

# Create a unit square Mesh. Let the mesh consist of 32 x 32 squares with each
# square divided into two triangles
mesh = UnitSquareMesh(32, 32)

# Define finite elements spaces and build mixed space
BDM = FiniteElement("BDM", mesh.ufl_cell(), 1)  # Order k for Stress
    # FunctionSpace; Brezzi-Douglas-Marini
DG  = FiniteElement("DG", mesh.ufl_cell(), 0)   # Order k - 1 for temperature
    # FunctionSpace; Discontinous Galerkin (Lagrange) elements
W = FunctionSpace(mesh, BDM * DG)
    # The second argument to FunctionSpace specifies underlying finite element,
    # where the mixed element is obtained by the * operator

# Specify the trial functions (the unknowns) and the test functions on W space
(sigma, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)

# Define the source function f. This is done just as for the Poisson demo:
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", \
                degree=2)

# Define the variational forms a and L. Since, u_0 = 0, the boundary term on
# the right-hand side (linear form) vanishes.
a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
L = - f*v*dx

# Boundary Conditions (BC)
#------------------------------------------------------------------------------
# Prescribe the boundary condition for the flux.

# Specifying the relevant part of the boundary can be done as for the Poisson
# demo (but now the top and bottom of the unit square is the essential boundary)

# Define essential boundary
def boundary(x):
    return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

# Define function G such that G \cdot n = g
class BoundarySource(UserExpression):
    def __init__(self, mesh, **kwargs):
        self.mesh = mesh
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        g = sin(5*x[0])
        values[0] = g*n[0]
        values[1] = g*n[1]
    def value_shape(self):
        return (2,)

G = BoundarySource(mesh, degree=2)

# Essential boundary conditions are specified through the class DirichletBC.
# Construct the essential BC by applying G to the first subspace of the mixed
# space (BDM).
bc = DirichletBC(W.sub(0), G, boundary)
    # Subspaces of a mixed FunctionSpace can be accessed by the method sub().

# Compute solution
w = Function(W)
solve(a == L, w, bc)
# The separate components of the solution can be extracted by using split()
(sigma, u) = w.split()

# Save solution in VTK format
file = XDMFFile("mixed_poisson.xdmf")
# Rename results for visualization in Paraview
sigma.rename("Stress", "sigma")
u.rename("Temperature", "u")
# Parameters will share the same mesh
file.parameters["flush_output"] = True
file.parameters["functions_share_mesh"] = True
# Save results
file.write(sigma, 0.0)
file.write(u, 0.0)

# # Plot the solutions (sigma and u) to examine the result.
# plt.figure()
# plot(sigma)
# plt.figure()
# plot(u)
# plt.show()
