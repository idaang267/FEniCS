# Import all modules
from dolfin import *

mesh = Mesh("2DShearTestRef.xml")
# subdomains = MeshFunction("size_t", mesh, "3d_Sphere_physical_region.xml")
# boundaries = MeshFunction("size_t", mesh, "3d_Sphere_facet_region.xml")

XDMFFile("mesh.xdmf").write(mesh)
# XDMFFile("mesh_physical_region.xdmf").write(subdomains)
# XDMFFile("mesh_facet_region.xdmf").write(boundaries)
