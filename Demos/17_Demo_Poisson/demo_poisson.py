# Poisson equation
# ===============
# This demo illustrates how to:
# * Solve a linear partial differential equation
# * Create and apply Dirichlet boundary conditions
# * Define Expressions
# * Define a Function Space
# * Create a SubDomain

# Implementation
#------------------------------------------------------------------------------
# First, import modules
from dolfin import *                # Module dolfin
import matplotlib.pyplot as plt     # Module matplotlib

# Create mesh of the domain and define finite element function space relative
# to mesh
mesh = UnitSquareMesh(32, 32) # Standard domain > use built-in mesh
    # class UnitSquareMesh: mesh consisting of 32 x 32 squares with each square
    # divided into two triangles
V = FunctionSpace(mesh, "Lagrange", 1)
    # FunctionSpace(mesh, finite element family, polynomial degree)
    # In this case, our space V consists of first-order, continuous Lagrange
    # finite element functions (continuous piecewise linear polynomials)

# Boundary Conditions
#------------------------------------------------------------------------------
# Consider the Dirichlet boundary condition. A simple Python function,
# returning a boolean, can be used to define the subdomain for the Dirichlet
# boundary condition (Gamma_D). The function should return True for those points
# inside the subdomain and False for the points outside. In our case, say that
# the points (x, y) such that x = 0 or x = 1 are on the inside of Gamma_D.

# Define Dirichlet boundary on the sides of the unit square (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Note that because of rounding-off errors, it is often wise to instead specify
# x < epsilon or x > 1 - epsilon where epsilon is a small number (such as
# machine precision)

# The Dirichlet boundary condition can be created using the class DirichletBC.
# DirichletBC(Function space the BC applies to,
#             value of BC,
#             part on the boundary on which the condition applies)

# Define boundary condition
u0 = Constant(0.0)
    # Value of the BC can be represented using a Constant, not an Expression
bc = DirichletBC(V, u0, boundary)
    # Function space is V and the Dirichlet boundary is defined above

# Express the variational problem. First, specify the trial function u and the
# test function v, both within the function space V. Define a TrialFunction and
# a TestFunction on the previously define FunctionSpace V.
u = TrialFunction(V)
v = TestFunction(V)

# The source f and the boundary normal derivative g are involved in the
# variational forms, and hence we must specify these. Both f and g are given by
# simple mathematical formulas, and can be easily declared using the Expression
# class.
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", \
                degree=2)
g = Expression("sin(5*x[0])", degree=2)

# Note that the strings defining f and g use C++ syntax since, for efficiency,
# DOLFIN will generate and compile C++ code for these expressions at run-time.

# Write down the bilinear form a and the linear form L (using UFL operators):
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Consider a solution of the variational problem. First, define a Function u
# within V, the finite element function space, to represent the solution.
u = Function(V)
    # Upon initialization, it is simply set to zero (u0)

# Next, call the solve function with the arguments a == L, u and bc as follows:
solve(a == L, u, bc)

# The function u will be modified during the call to solve. The default settings
# for solving a variational problem have been used. However, the solution
# process can be controlled in much more detail if desired.

# A Function can be manipulated in various ways, in particular, it can be
# plotted and saved to file. Here, output the solution to a VTK file (using the
# suffix .pvd) for later visualization and also plot it using the plot command:

# Save solution in VTK format
file = File("poisson.pvd")
# Rename results for visualization in Paraview
u.rename("Displacement", "u")
file << u

# Plot solution
plot(u)
plt.show()
