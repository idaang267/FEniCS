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
# Uses the same definitions of functions and boundaries as in the demo for
# Poisson's equation.

# First, the required modules are imported
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

# Define boundary classes for Dirichlet BC
class BotBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS

class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 1.0 - DOLFIN_EPS

BotBoundary = BotBoundary()
TopBoundary = TopBoundary()

# Number of elements across
N = 20
# For saving
simulation_params = "Results_N_%.0f" % (N)
savedir = "../Results/PoissonMixed/" + simulation_params + "/"

# Create a unit square Mesh, N x N, where each square is divided into two triangles
mesh = UnitSquareMesh(N, N)

SubDomains = MeshFunction("size_t", mesh, mesh.topology().dim())
Lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
Points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

# show lines of interest
Lines.set_all(0)
BotBoundary.mark(Lines,1)
TopBoundary.mark(Lines,2)
file_results = XDMFFile(savedir + "/" + "Lines.xdmf")
file_results.write(Lines)

#makes sure each subdomain has values specified by mark
ds = Measure('ds', domain=mesh, subdomain_data=Lines)

# Define finite elements spaces for mixed function space definition
#------------------------------------------------------------------------------
# FunctionSpace; Brezzi-Douglas-Marini, Order k for stress
BDM = FiniteElement("BDM", mesh.ufl_cell(), 1)
# FunctionSpace; Discontinous Galerkin (DG), Order k - 1 for temperature
DG  = FiniteElement("DG", mesh.ufl_cell(), 0)
# The mixed element is obtained by the * operator
W = FunctionSpace(mesh, BDM * DG)

# Specify the trial functions (the unknowns) and the test functions on W space
(q, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)
w = Function(W)                 # Solutions for (heat flux, temperature)

# Boundary Conditions (BC) and expressions
#------------------------------------------------------------------------------
# Define the expressions for the Neumann BC which now become Dirichlet BCs
g0 = Expression(("0", "sin(5*x[0])"), degree=2)
g1 = Expression(("0", "-sin(5*x[0])"), degree=2)
# Prescribe the boundary condition on the boundary, top and bottom
BotBC = DirichletBC(W.sub(0), g0, BotBoundary)
TopBC = DirichletBC(W.sub(0), g1, TopBoundary)

# bc1 = DirichletBC(W.sub(0), g1, Boundary)
bc = [BotBC, TopBC]

# Define the source function, which is the same as the Poisson demo:
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", \
                degree=2)
# Define the variational forms a and L
a = (dot(q, tau) + div(tau)*u + div(q)*v)*dx
L = - f*v*dx
# Solve
solve(a == L, w, bc)

# The separate components of the solution can be extracted by using split()
(q, u) = w.split()

# Save solution in XDMF format
file = XDMFFile(savedir + "/Poisson.xdmf")
# Parameters will share the same mesh
file.parameters["flush_output"] = True
file.parameters["functions_share_mesh"] = True

# Rename results for visualization in Paraview
q.rename("Heat Flux", "q")
u.rename("Temperature", "u")

# Save results
file.write(q, 0.0)
file.write(u, 0.0)

# Plot the solutions (q and u) to examine the result.
plt.figure()
plot(q)
plt.savefig(savedir + "/HeatFlux.pdf", transparent=True, bbox_inches='tight')
plt.figure()
plot(u)
plt.savefig(savedir + "/Temperature.pdf", transparent=True, bbox_inches='tight')
plt.show()
