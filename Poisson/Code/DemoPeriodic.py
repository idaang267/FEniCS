# This demo program solves Poisson's equation
#
#     - div grad u(x, y) = f(x, y)
#
# on the unit square with homogeneous Dirichlet boundary conditions at y = 0, 1
# and periodic boundary conditions at x = 0, 1.
#
# Original implementation by Anders Logg
# Copyright (C) 2007 Kristian B. Oelgaard
#
# Modified by Ida Ang in 2023
#
# First added:  2007-11-15
# Last changed: 2023-04-03

# Import Modules
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Source term
class Source(UserExpression):
    def eval(self, values, x):
        dx = x[0] - 0.5
        dy = x[1] - 0.5
        values[0] = x[0]*sin(5.0*DOLFIN_PI*x[1]) + 1.0*exp(-(dx*dx + dy*dy)/0.02)
    def value_shape(self):
         return ()

# Subdomain for Dirichlet boundary condition, top and bottom
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < 0.0 + DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

# Subdomain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]

TopBot = DirichletBoundary()
PBC = PeriodicBoundary()

# For saving
savedir = "../Results/PoissonPeriodic/"

# Create mesh and visualize boundaries
N = 32
mesh = UnitSquareMesh(N, N)
SubDomains = MeshFunction("size_t", mesh, mesh.topology().dim())
Lines = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# show lines of interest
Lines.set_all(0)
TopBot.mark(Lines,1)
PBC.mark(Lines,2)
file_results = XDMFFile(savedir + "Lines.xdmf")
file_results.write(Lines)

# Set up function spaces and define variational problem
V = FunctionSpace(mesh, 'CG', 2, constrained_domain=PBC)
DG0 = FunctionSpace(mesh,'DG',0)

# Trial, test, and solution definitions
u = TrialFunction(V)
v = TestFunction(V)
Sol = Function(V)

# Create Dirichlet boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, TopBot)
f = Source(degree=2)

# Define bilinear and linear form
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
solve(a == L, Sol, bc)

# Save solution to file
file = XDMFFile(savedir + f"./results.xdmf")
# Saves the file in case of interruption
file.parameters["rewrite_function_mesh"] = False
file.parameters["functions_share_mesh"]  = True
# Project source term to discontinuous galerkin space space
SolSource = project(f, DG0)
# Rename results for visualization in Paraview
Sol.rename("Potential Field", "Sol")
file.write(Sol, 0.0)
file.write(SolSource, 0.0)

# Postprocessing
Points = np.linspace(0, 1.0, N)         # Points along the profile
SolLeft = np.zeros((1, len(Points)))    # Array for storing sol along profile
SolRight = np.zeros((1, len(Points)))
SolLeftIn = np.zeros((1, len(Points)))
SolRightIn = np.zeros((1, len(Points)))

for (i, int) in enumerate(Points):
    SolLeft[:, i] = Sol(0,int)
    SolLeftIn[:, i] = Sol(0.1,int)
    SolRight[:, i] = Sol(1,int)
    SolRightIn[:, i] = Sol(0.9,int)

# Plot and save solution to figure
font = {'family': 'serif', 'weight': 'normal', 'size': 16}
FontLegend = {'family': 'serif', 'weight': 'normal', 'size': 14}
plt.rc('font', **font)

# Figure 1
fig, ax = plt.subplots(1,2, figsize=(11.5,5), sharey=True)
ax[0].plot(Points, SolLeft[0,:],"r", linewidth=2, label="Left")
ax[0].plot(Points, SolRight[0,:], "k--", linewidth=2, label="Right")
ax[1].plot(Points, SolLeftIn[0,:],"r", linewidth=2, label="0.1 from Left")
ax[1].plot(Points, SolRightIn[0,:], "k--", linewidth=2, label="0.1 from Right")

for i in [0,1]:
    ax[i].set_xlabel("Y Profile", fontdict=font)
    ax[i].set_ylabel("Potential Field", fontdict=font)
    ax[i].set_xlim([0, 1])
    ax[i].set_ylim([0, 0.015])
    ax[i].legend(loc="best", frameon=True, prop=FontLegend)

plt.savefig("../Images/Periodic/PeriodicSolProfile.pdf", transparent=True, bbox_inches='tight')
plt.close()

# Figure 2
plt.figure()
plot(SolSource)
plt.title("Source Function")
plt.xlabel("X Coordinate", fontdict=font)
plt.ylabel("Y Coordinate", fontdict=font)
plt.savefig("../Images/Periodic/PeriodicSource.pdf", transparent=True, bbox_inches='tight')
plt.close()

# Figure 3 
plt.figure()
plot(Sol)
plt.title("Temperature")
plt.xlabel("X Coordinate", fontdict=font)
plt.ylabel("Y Coordinate", fontdict=font)
plt.savefig("../Images/Periodic/PeriodicTemp.pdf", transparent=True, bbox_inches='tight')
plt.close()
