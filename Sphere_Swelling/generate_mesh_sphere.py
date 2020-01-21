# Copyright (C) 2016-2019 by the multiphenics authors
#
# This file is part of multiphenics.
#
# multiphenics is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiphenics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiphenics. If not, see <http://www.gnu.org/licenses/>.
#

from dolfin import *
from mshr import *
from multiphenics import *

# Create mesh
domain = Sphere(Point(0.0, 0.0, 0.0), 1.0)
mesh = generate_mesh(domain, 15)

# Create subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), mesh.domains())

# Create boundaries
class OnBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
class OnInterface(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
on_boundary = OnBoundary()
on_boundary.mark(boundaries, 1)
on_interface = OnInterface()
on_interface.mark(boundaries, 2)

# Create restriction
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[1], l0) and on_boundary)

boundary_restriction = MeshRestriction(mesh, on_boundary)
top = Top()
top_restriction = MeshRestriction(mesh, top)
interface_restriction = MeshRestriction(mesh, on_interface)

# Save
File("sphere.xml") << mesh
File("sphere_physical_region.xml") << subdomains
File("sphere_facet_region.xml") << boundaries
XDMFFile("sphere.xdmf").write(mesh)
XDMFFile("sphere_physical_region.xdmf").write(subdomains)
XDMFFile("sphere_facet_region.xdmf").write(boundaries)

File("sphere_restriction_boundary.rtc.xml") << boundary_restriction
File("sphere_restriction_interface.rtc.xml") << interface_restriction
XDMFFile("sphere_restriction_boundary.rtc.xdmf").write(boundary_restriction)
XDMFFile("sphere_restriction_interface.rtc.xdmf").write(interface_restriction)
