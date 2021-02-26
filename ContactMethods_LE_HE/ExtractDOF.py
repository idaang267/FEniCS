# Import modules
from dolfin import *
import numpy as np

# Square unit mesh
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), 10, 10, 10)
# Lagrange function space of degree 1
V = FunctionSpace(mesh, 'CG', 1)
# Test expression
exp = Expression('x[0]*x[1]', degree=2)
# Interpolate test expression to function space
u = interpolate(exp, V)
# Returns topological dimension of function space, V
dim = V.dim()   # Calculated as (# + 1) x (# + 1)
# Gives dimensions: 2D vs 3D
N = mesh.geometry().dim()
#
dofmap = V.dofmap()
dof = V.tabulate_dof_coordinates().reshape(dim, N)
x = dof[:, 0]
y = dof[:, 1]
# Extract dof indices where some condition is met
#indices = np.where(np.logical_and(x > 0.5, y == 1.0))[0]
indices = np.where(np.logical_and(x > 0.45, x < 0.55))[0]
# Get coordinates of dof
xs = dof[indices]
# Get value of dof
vals = u.vector()[indices]

for x, v in zip(xs, vals):
    print(x, v)
    if x[1] == 0.5 and x[2] == 0.5:
        print(v)
