# Copyright (C) 2016-2019 by the multiphenics authors
#
# This file is part of multiphenics.

'''
In this example, a Laplace problem is solved with non-homogeneous Dirichlet
boundary conditions (BCs), which are imposed with Lagrange multipliers. Standard
FEniCS code does not easily support Lagrange multipliers, because FEniCS does
not support subdomain/boundary restricted function spaces. One would have to
declare the Lagrange multiplier on the entire domain and constrain it in the
interior. This procedure would require the definition of suitable MeshFunctions
to constrain the additional DOFs, resulting in a (1) cumbersome mesh definition
for the user and (2) unnecessarily large linear system. This task is more easily
handled by multiphenics by providing a restriction in the definition of the
(block) function space. Such a restriction (which is basically a collection of
MeshFunctions) can be generated from a SubDomain object.
'''

# Import all modules
from dolfin import *
from multiphenics import *
from mshr import *
from ufl import rank

snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "mumps",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": False}}

# Parameters
depth = 0.1
# Steps
steps = 0
tot_steps = 5
# Define the indenter (obstacle) with a parabola
R = 0.25                            # Radius
depth = 0.01                        # Depth

# Create Mesh
domain = Rectangle(Point(0., 0.), Point(1.0, 1.0))
mesh = generate_mesh(domain, 15)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())
# Create boundaries
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Define class
class top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and on_boundary

class bot(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

# Mark boundary to restrict the lagrange multiplier to
top = top()
bot = bot()
top.mark(boundaries, 1)

# Dirichlet boundary
boundary_restriction = MeshRestriction(mesh, top)

# Function space
V = VectorFunctionSpace(mesh, "Lagrange", 2)
# Block mixed equal order function space
W = BlockFunctionSpace([V, V], restrict=[None, boundary_restriction])

# Define trial and test functions
du = BlockTrialFunction(W)
ul = BlockFunction(W)
v = BlockTestFunction(W)
# Split trial and test functions into displacement and lagrange multiplier
(u, l) = block_split(ul)
(v_u, v_l) = block_split(v)

# Measure
dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

indent = Expression(("0.0","0.0"), depth=depth, degree=1)
bc_bot = DirichletBC(W.sub(0), indent, bot)
bc = BlockDirichletBC([bc_bot, None])

# Assemble
# Type of source term?
f = Expression(("0.0","0.0"), element=V.ufl_element())
g = Expression(("0.0","-depth+(pow(x[0],2))/(2*R)"), depth=depth, R=R, element=V.ufl_element())
F = [inner(grad(u), grad(v_u))*dx + inner(l,v_u)*ds + inner(f,v_u)*dx,
    dot(v_l,(u-g))*ds]

J = block_derivative(F, ul, du)
problem = BlockNonlinearProblem(F, ul, bc, J)
solver = BlockPETScSNESSolver(problem)
solver.parameters.update(snes_solver_parameters["snes_solver"])

results = XDMFFile("lagrangeLinear.xdmf")
while steps <= tot_steps:
    print(steps)
    solver.solve()
    # Update
    steps += 1
    depth += 0.01
    g.depth = depth
    # Split results
    (u, l) = ul.block_split()
    # Save results
    u.rename("Displacement", "u")
    l.rename("Lagrange Multiplier", "l")
    # Parameters will share the same mesh
    results.parameters["flush_output"] = True
    results.parameters["functions_share_mesh"] = True
    # Write
    results.write(u, steps)
    results.write(l, steps)
