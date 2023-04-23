# Heat equation with Dirichlet conditions.

# Test problem is chosen to give an exact solution at all nodes of the mesh
#  u = 1 + x^2 + alpha*y^2 + beta*t
#  f = beta - 2 - 2*alpha

from __future__ import print_function
from fenics import *
import numpy as np

# Constants
#------------------------------------------------------------------------------
T = 10.0            # Final time
NumSteps = 50      # Number of time steps
dt = T / NumSteps  # Time step size
# Constants used for contruction of the test problem
alpha = 3          # Parameter alpha
beta = 1.2         # Parameter beta
nx = ny = 40        # Mesh resolution

# Create mesh and define function space
#------------------------------------------------------------------------------
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'P', 1)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)    # Found from the guess for u
SimPar = "S_%.0d_N_%.0d" % (NumSteps, nx)
Dir = "../Results/Interpolate/" + SimPar

# Define initial condition, (IC), value at the previous time step
#------------------------------------------------------------------------------
# Can be computed by either projection or interpolation of u_o
# NOTE: projection results in approximate values while interpolation is exact

# Expression given: u_D = u = 1 + x^2 + alpha*y^2 + beta*t
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t',
                 degree=2, alpha=alpha, beta=beta, t=0)

# Since we set t=0 for the boundary value u_D, we can use u_D to specify
u_n = interpolate(u_D, V)
# u_n = project(u_D, V)

# Define boundary conditions (BC)
#------------------------------------------------------------------------------
# Essential boundary conditions along the entire boundary
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Weak form where we define F
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
# FEniCS can determine which is the bilinear form, a, and linear form, L.
a, L = lhs(F), rhs(F)

# Save results to an .xdmf file since we have multiple fields, result is
# time dependent (displacement, time)
file_results = XDMFFile(Dir + "/Results.xdmf")
# Parameters will share the same mesh
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

# Time-stepping loop
u = Function(V)         # We set this as a trial function as well as a function

# Postprocessing
PostTxt = np.zeros((NumSteps, 3))

t = 0                   # Define time initialization
for n in range(NumSteps):
    # Compute solution
    solve(a == L, u, bc)

    u_e = interpolate(u_D, V)
    u.rename("Temperature", "u")
    u_e.rename("Temp Exact", "u_e")

    # Write solution
    file_results.write(u,t)
    file_results.write(u_e,t)

    # Compute error at vertices
    error = np.abs(u_e.vector() - u.vector()).max()
    print('t = %.2f: error = %.3g' % (t, error))

    PostTxt[n, :] = np.array([n, t, error])
    np.savetxt(Dir + '/PostProc.txt', PostTxt)

    # Update current time
    t += dt
    u_D.t = t           # Update time within expression class
    u_n.assign(u)       # Update previous solution
