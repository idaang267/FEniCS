# Copyright (C) 2016-2019 by the multiphenics authors
#
# This file is part of multiphenics.
#
# multiphenics is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# multiphenics is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiphenics. If not, see <http://www.gnu.org/licenses/>.

"""
In this tutorial we compare the formulation and solution of a Navier-Stokes by
standard FEniCS code (using the MixedElement class) and multiphenics code.
"""

# Solver parameters
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": True}}

# Import Modules
from dolfin import *
from multiphenics import *
from mshr import *
# Import specific functions from modules
from numpy import isclose
from ufl import rank

# Parameters
steps = 0                           # Steps (update within loop)
tot_steps = 5                       # Total steps
# Constitutive Parameters
dyn_vis = Constant(0.01)            # Dynamic Viscosity

# Create mesh
domain = Box(Point(0.,0.,0), Point(1.0,1.0,1.0))
mesh = generate_mesh(domain, 10)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, 2)
subdomains.set_all(0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
top = Top()
bot = Bottom()
top.mark(boundaries, 1)

boundary_restriction = MeshRestriction(mesh, top)

# Measure
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

# Function spaces
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange",2)
W = BlockFunctionSpace([V, Q, Q], restrict=[None, None, boundary_restriction])

# Test and trial functions
dup = BlockTrialFunction(W)
v = BlockTestFunction(W)
up = BlockFunction(W)
# Split
(v_u, v_p, v_l) = block_split(v)
(u, p, l) = block_split(up)

# Constitutive parameters for boundary conditions
fix = Constant((0., 0., 0.))

# Boundary conditions
fixed_bot = DirichletBC(W.sub(0), fix, bot)
bc = BlockDirichletBC([[fixed_bot], []])

# Weak forms
depth = 0.00
u_in = Expression(("-depth"), depth = depth, degree = 0)
n   = FacetNormal(mesh)

F = [dyn_vis*inner(grad(u), grad(v_u))*dx - p*div(v_u)*dx + l*v_u[1]*ds,
     v_p*div(u)*dx,
     v_l*(u[1]-u_in)*ds]
J = block_derivative(F, up, dup)

# Solve
problem = BlockNonlinearProblem(F, up, bc, J)
solver = BlockPETScSNESSolver(problem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

# Write to file
results = XDMFFile("b_navier_stokes.xdmf")
while steps <= tot_steps:
    solver.solve()
    # Update
    steps += 1
    depth += 0.01
    u_in.depth = depth
    # Extract solutions
    (u, p, l) = up.block_split()
    # Rename
    u.rename("Displacement", "u")
    p.rename("Pressure", "p")
    l.rename("Lagrange Multiplier","l")
    # Parameters will share the same mesh
    results.parameters["flush_output"] = True
    results.parameters["functions_share_mesh"] = True
    # Write to .xdmf results file
    results.write(u, steps)
    results.write(p, steps)
    results.write(l, steps)
