'''
Enforcing Dirichlet boundary conditions weakly using Lagrange multipler techni-
que introduced in
    [Babuska; The finite element method with lagrangian multipliers]

    Simple Poisson problem:
        -laplace(u) = f in [0, 1]^2
                  u = g on boundary
'''

from collections import namedtuple
from dolfin import *

Result = namedtuple('Result', ['h', 'H1', 'L2'])

u_exact = Expression('x[0]*x[1] + sin(pi*x[0])*cos(2*pi*x[1])', degree=2)
f = Expression('5*pi*pi*sin(pi*x[0])*cos(2*pi*x[1])', degree=2)

# Definition for solving the Poisson problem
def lagrange_solver(N):

    mesh = UnitSquareMesh(N, N)

    # Equal order mixed function space for displacement and lagrange multiplier
    V = FunctionSpace(mesh, 'CG', 1)
    Velem = V.ufl_element()
    W = FunctionSpace(mesh, Velem*Velem)

    # Define functions in mixed function space W
    u, lmbda = TrialFunctions(W)    # Displacement and lagrange multiplier
    v, mu = TestFunctions(W)

    # Weak form
    a = inner(grad(u), grad(v))*dx + inner(lmbda, v)*ds + inner(mu, u)*ds
    L = inner(f, v)*dx + inner(mu, u_exact)*ds

    A, b = assemble_system(a, L)
    # Eliminate dofs of W.sub(1) which are unset since the form has ds, not dx.
    A.ident_zeros()

    Uh = Function(W)
    solve(A, Uh.vector(), b)

    uh, lmbdah = Uh.split(True)
    # plot(uh, interactive=True)
    file_results = XDMFFile("lagrange_multiplier.xdmf")
    # Rename results for visualization in Paraview
    uh.rename("Displacement", "uh")
    lmbdah.rename("Lagrange Multiplier", "lmbdah")
    # Parameters will share the same mesh
    file_results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    file_results.write(uh, 0.0)
    file_results.write(lmbdah, 0.0)

    # Error norms
    error_L2 = errornorm(u_exact, uh, 'L2')
    error_H1 = errornorm(u_exact, uh, 'H1')
    return Result(h=mesh.hmin(), H1=error_H1, L2=error_L2)

# -----------------------------------------------------------------------------

norm_type = 'L2'
# Input into lagrange_solver the amount of elements
R = lagrange_solver(N=4)
h_ = R.h
e_ = getattr(R, norm_type)
for N in [8, 16, 32, 64, 128]:
    R = lagrange_solver(N=N)
    h = R.h
    e = getattr(R, norm_type)
    rate = ln(e/e_)/ln(h/h_)
    print('{h:.3E} {e:.3E} {rate:.2f}'.format(h=h, e=e, rate=rate))
