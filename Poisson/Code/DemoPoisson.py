# Poisson equation
# ===============
# Implementation
# First, import modules
from dolfin import *                # Module dolfin
import matplotlib.pyplot as plt     # Module matplotlib
import numpy as np

# Define boundary classes for Dirichlet BC which specifies left and right side
# of the square
class boundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < -0.5 + DOLFIN_EPS or x[0] > 0.5 - DOLFIN_EPS

# All parameters
# -----------------------------------------------------------------------------
N = 20            # Number of elements across

# The source f definition
f = Expression("10*exp(-(pow(x[0], 2) + pow(x[1], 2)) / 0.02)", degree=1)
# Boundary normal derivative g
g_x = Expression("sin(5*x[0])", degree=2)

# For saving
simulation_params = "Results_N_%.0f" % (N)
savedir = "../Results/Poisson/" + simulation_params + "/"

# Create mesh of the domain and define finite element function space
mesh = RectangleMesh(Point(-0.5, -0.5), Point(0.5, 0.5), int(N), int(N))

SubDomain = MeshFunction("size_t", mesh, mesh.topology().dim())
Lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
points = MeshFunction("size_t", mesh, mesh.topology().dim() - 2)

# show lines of interest
Lines.set_all(0)
boundary().mark(Lines,1)
file = XDMFFile(savedir + "/" + "Lines.xdmf")
file.write(Lines)

#makes sure each subdomain has values specified by mark
ds = Measure('ds', domain=mesh, subdomain_data=Lines)

# Function spaces
V = FunctionSpace(mesh, "Lagrange", 1)
u = TrialFunction(V)
v = TestFunction(V)
# Define a Function within V, the FE function space, to represent the solution.
Sol = Function(V)

# Boundary Conditions
#------------------------------------------------------------------------------
# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary())

# Write down the bilinear form a and the linear form L (using UFL operators):
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g_x*v*ds

# Next, call the solve function with the arguments a == L, u and bc as follows:
solve(a == L, Sol, bc)

# Save results to an .xdmf file since we have multiple fields
file = XDMFFile(savedir + f"./results.xdmf")
file.parameters["rewrite_function_mesh"] = False
file.parameters["flush_output"] = True

# Rename results for visualization in Paraview
Sol.rename("Potential Field", "Sol")
file.write(Sol)

coords = V.tabulate_dof_coordinates()
vec = Sol.vector().get_local()
PostProc = np.zeros((len(coords), 3))
Count = 0
for Coord, Val in zip(coords, vec):
    PostProc[Count] = np.array([Coord[0], Coord[1], Val])
    Count += 1
np.savetxt(savedir + f'./results.txt', PostProc)

# Plot the solutions (q and u) to examine the result.
plt.figure()
plot(Sol)
plt.savefig("../Images/Poisson/Temperature.pdf", transparent=True, bbox_inches='tight')
plt.close()
