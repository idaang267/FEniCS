'''
Enforcing Dirichlet boundary conditions via penalty term following

    [Babuska; The finite element method with lagrangian multipliers]

    Simple Poisson problem:
        -laplace(u) = f in [0, 1]^2
                  u = g on boundary
'''

from collections import namedtuple
from dolfin import *
from ufl import *

Result = namedtuple('Result', ['h', 'H1', 'L2'])

u_exact = Expression('x[0]*x[1] + sin(pi*x[0])*cos(2*pi*x[1])', degree=2)
f = Expression('5*pi*pi*sin(pi*x[0])*cos(2*pi*x[1])', degree=2)

# Definition for solving the Poisson problem
def penalty_solver(N, penalty):
    # Define square mesh
    mesh = UnitSquareMesh(N, N)

    # Function space for displacement
    V = FunctionSpace(mesh, 'CG', 1)
    # Define trial and test function
    u = TrialFunction(V)
    v = TestFunction(V)

    h = mesh.ufl_cell().MaxFacetEdgeLength

    # Penalty Parameter
    sigma = Constant(penalty)

    # Weak Form
    a = inner(grad(u), grad(v))*dx + (1/h**sigma)*inner(u, v)*ds
    L = inner(f, v)*dx + (1/h**sigma)*inner(u_exact, v)*ds

    A, b = assemble_system(a, L)

    uh = Function(V)

    # Write results
    file_results = XDMFFile("penalty.xdmf")
    uh.rename("Displacement", "uh")
    file_results.write(uh, 0.0)

    # Solve
    solve(A, uh.vector(), b)

    # Error norms
    error_L2 = errornorm(u_exact, uh, 'L2')
    error_H1 = errornorm(u_exact, uh, 'H1')

    return Result(h=mesh.hmin(), H1=error_H1, L2=error_L2)

# -----------------------------------------------------------------------------

norm_type = 'H1'
penalty_value = 2  #[1..10] all okay, f has nice regularity
# Input into lagrange_solver the amount of elements and penalty parameter value
R = penalty_solver(N=4, penalty=penalty_value)
h_ = R.h
e_ = getattr(R, norm_type)
for N in [8, 16, 32, 64, 128]:
    R = penalty_solver(N=N, penalty=penalty_value)
    h = R.h
    e = getattr(R, norm_type)
    rate = ln(e/e_)/ln(h/h_)
    print('{h:.3E} {e:.3E} {rate:.2f}'.format(h=h, e=e, rate=rate))
