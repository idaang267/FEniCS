from dolfin import *

# Define mesh and subdomain:
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), 30, 30, 20)
d = mesh.topology().dim()
#filmx = CompiledSubDomain("x[1] > 1.0-x[0]-DOLFIN_EPS")
#domain = CompiledSubDomain("pow(x[0]-0.5,2) <= pow(0.3,2) - pow(x[2]-0.5,2) - DOLFIN_EPS" or "x[1] > 0.5")

# Define class

refDomain = refDomain()
# Refinement using the `refine()` function and a Boolean `MeshFunction`:
r_markers = MeshFunction("bool", mesh, d, False)
refDomain.mark(r_markers, True)
refinedMesh = refine(mesh,r_markers)

# Transfering a non-negative integer-valued (`size_t`) `MeshFunction` to the
# refined mesh using the `adapt()` function:
meshFunctionToAdapt = MeshFunction("size_t", mesh, d, 0)
refDomain.mark(meshFunctionToAdapt,1)
adaptedMeshFunction = adapt(meshFunctionToAdapt,refinedMesh)

# Plot results:
File("Mesh.pvd") << adaptedMeshFunction
