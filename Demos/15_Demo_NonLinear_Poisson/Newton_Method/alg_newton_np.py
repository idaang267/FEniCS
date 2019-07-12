# Nonlinear Poisson Equation with Newton Algorithm
# ==========================
# Dirichlet conditions in x-direction and homogeneous Neumann (symmetry)
# conditions in all other directions. The domain is the unit hypercube in a
# given dimension.
#
# -div(q(u)*nabla_grad(u)) = 0 where q(u) = (1+u)^m
# u = 0 at x=0, u=1 at x=1,     on gamma_D
# du/dn=0                       on gamma_N
# NOTE: q(u) is different from the other non-linear example

# Solution method: Newton

# First, modules are imported
from dolfin import *                # Module dolfin
import numpy, sys
import matplotlib.pyplot as plt     # Module matplotlib

# Create mesh and define function space
degree = int(sys.argv[1])
divisions = [int(arg) for arg in sys.argv[2:]]
d = len(divisions)
domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
mesh = domain_type[d-1](*divisions)
V = FunctionSpace(mesh, 'Lagrange', degree)

# Define boundary conditions for initial guess
tol = 1E-14
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol
def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol

Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary) # u = 0 at x = 0
Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary) # u = 1 at x = 1
bcs = [Gamma_0, Gamma_1]

# Define variational problem for initial guess (q(u)=1, i.e., m=0)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(nabla_grad(u), nabla_grad(v))*dx
f = Constant(0.0)
L = f*v*dx
A, b = assemble_system(a, L, bcs)
u_k = Function(V)
solve(A, u_k.vector(), b, 'lu')

# All Dirichlet conditions must be zero for the correction function in a Newton
# method: u = u_k + omega*du. In order for u^k to fulfill the Dirichlet
# conditions for u, du most be zero on the boundaries
Gamma_0_du = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1_du = DirichletBC(V, Constant(0.0), right_boundary)
bcs_du = [Gamma_0_du, Gamma_1_du]

m = 2 # Choice of nonlinear coefficient
def q(u):
    return (1+u)**m
def Dq(u):
    return m*(1+u)**(m-1)

# Define variational problem for the matrix and vector in a Newton iteration
du = TrialFunction(V) # du is delta u where u = u_k + omega*du
a = inner(q(u_k)*nabla_grad(du), nabla_grad(v))*dx + \
    inner(Dq(u_k)*du*nabla_grad(u_k), nabla_grad(v))*dx
L = -inner(q(u_k)*nabla_grad(u_k), nabla_grad(v))*dx

# Newton iteration at the algebraic level
du = Function(V)
u  = Function(V)
omega = 1.0         # relaxation parameter
eps = 1.0           # error measure
tol = 1.0E-5        # tolerance
iter = 0            # iteration counter
maxiter = 25        # Max # of iterations allowed
while eps > tol and iter < maxiter:
    iter += 1
    print(iter, 'iteration')
    A, b = assemble_system(a, L, bcs_du) # u_k must have the correct BCs
    solve(A, du.vector(), b)
    eps = numpy.linalg.norm(du.vector().array(), ord=numpy.Inf)
    print('Norm:', eps)
    u.vector()[:] = u_k.vector() + omega*du.vector()
    # or
    #u.vector()[:] += omega*du.vector()
    # or
    #u.assign(u_k)  # u = u_k
    #u.vector().axpy(omega, du.vector())
    u_k.assign(u)

# Print out statements
convergence = "convergence after %d Newton iterations at the algebraic level" % iter
if iter >= maxiter:
    convergence = 'no ' + convergence

print('Solution of the nonlinear Poisson problem div(q(u)*nabla_grad(u)) = f, \
with f=0, q(u) = (1+u)^m, u=0 at x=0 and u=1 at x=1. %s %s' \
% (mesh,convergence))

# Find max error
# u_exact = Expression("pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1", m=m)
# u_e = interpolate(u_exact, V)
# diff = numpy.abs(u_e.vector().array() - u.vector().array()).max()
# print('Max error:', diff)

# Plot Solution
plt.figure()
plot(u, title="Solution") # Plot solution
plt.show()
