# Nonlinear Poisson equation
# ==========================

# This demo illustrates how to solve a nonlinear partial differential equation
# (in this case a nonlinear variant of Poisson's equation)

# First, the dolfin and matplotlib modules are imported
from dolfin import *
import matplotlib.pyplot as plt

# Boundary Conditions
# --------------
# A simple Python function, returning a boolean, can be used to define the
# subdomain for the Dirichlet boundary condition (BC) Gamma_D. The function
# should return True for those points inside the subdomain and False for the
# points outside. In our case, we want to say that the points (x, y) such that
# x = 1 are on the inside of Gamma_D.

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary
    # Note: because of rounding-off errors, it is often wise to instead specify
    # |x - 1| < epsilon, where epsilon is a small number
    # DOLFIN_EPS: machine precision?

# We then define a mesh of the domain and a finite element function space V
# relative to this mesh.

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32) # Built-in mesh
    # Consists of 32 x 32 squares with each square divided into two triangles
V = FunctionSpace(mesh, "CG", 1)
    # FunctionSpace(mesh, finite element family, polynomial degree)
    # In this case, we use 'CG', for Continuous Galerkin, as a synonym for
    # 'Lagrange'. With degree 1, we get the standard linear Lagrange element
    # otherwise known as continuous piecewise linear polynomials, which are
    # triangles with nodes at the three vertices.

# The Dirichlet boundary condition can be created using the class DirichletBC.
# DirichletBC(Function space the BC applies to,
#             value of BC,
#             part on the boundary on which the condition applies)

# Define boundary condition on Gamma_D (Dirichlet Boundary)
u_0 = Constant(1.0)
    # Value of the BC can be represented using a Constant, not an Expression
bc = DirichletBC(V, u_0, DirichletBoundary())
    # Function space is V and the Dirichlet boundary is defined above

# Next, we want to express the variational problem. First, we need to specify
# the function u which represents the solution. Upon initialization, it is
# simply set to the zero function, which will represent the initial guess u_0.
# We do this by defining a Function and a TestFunction on the previously
# defined FunctionSpace V.

# A Function represents a function u living in a finite element function space.
u = Function(V)
# The Test Function v is specified, also living in the function space V.
v = TestFunction(V)

# Further, the source f is involved in the variational forms, and must be
# specified. We have f given by a simple mathematical formula, which can be
# easily declared using the Expression class.
f = Expression("x[0]*sin(x[1])", degree=2)

# By defining the function in this step and omitting the trial function we tell
# FEniCS that the problem is nonlinear. With these ingredients, we can write
# down the semilinear form F. Use Unifed Form Language (UFL) operators
q = 1 + u**2
F = inner(q*grad(u), grad(v))*dx - f*v*dx

# Now, we have specified the variational forms and can consider the solution of
# the variational problem. Next, we can call the solve function with the
# arguments F == 0, u, bc and solver parameters as follows

# Compute solution
solve(F == 0, u, bc,
      solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})

# The Newton procedure is considered to have converged when the residual r_n at
# iteration n is less than the absolute tolerance or the relative residual
# \frac{r_n}{r_0} is less than the relative tolerance.

# A Function can be manipulated in various ways, in particular, it can be
# plotted and saved to file.

plt.figure()
plot(u, title="Solution") # Plot solution
plt.figure()
plot(grad(u), title="Solution gradient") # Plot solution gradient
plt.show()

# Output the solution to a VTK file (using the suffix .pvd) for later
# visualization
File("nonlinear_poisson_mesh.pvd") << mesh
File("nonlinear_poisson.pvd") << u
