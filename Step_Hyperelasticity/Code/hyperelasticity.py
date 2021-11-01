# Hyperelasticity
# ===============
# This example demonstrates the solution of a three-dimensional elasticity
# problem.

# First, the required modules are imported
from dolfin import *
import matplotlib.pyplot as plt

# The behavior of the FEniCS Form Compiler, FFC, can be adjusted by prescribing
# various parameters. Use the UFL Analyser and Compiler System, UFLACS, backend
# of FFC.

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
    # UFLACS is a collection of algorithms for processing symbolic unified form
    # languages forms and expressions

# Create tetrahedral mesh of the domain and a function space on this mesh
mesh = UnitCubeMesh(24, 16, 16) # Domain - Omega
    # A unit cube mesh with 25 ( = 24 + 1) vertices in one direction and
    # 17 ( = 16 + 1) vertices in the other two direction.
V = VectorFunctionSpace(mesh, "Lagrange", 1)
    # On this mesh, we define a function space of continuous piecewise linear
    # vector polynomials (a Lagrange vector element space)

    # Note that VectorFunctionSpace creates a function space of vector fields.
    # The dimension of the vector field (the number of components) is assumed
    # to be the same as the spatial dimension, unless otherwise specified.

# The portions of the boundary on which Dirichlet boundary conditions will be
# applied are defined

# Mark boundary subdomains
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
    # Gamma_{D_0}: The boundary subdomain 'left' corresponds to the part of the
    # boundary on which 'x=0'
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
    # Gamma_{D_1}: The boundary subdomain 'right' corresponds to the part of the
    # boundary on which 'x=1'

    # Note: C++ syntax is used in the `CompiledSubDomain` function since the
    # function will be automatically compiled into C++ code for efficiency.

    # The (built-in) variable `on_boundary` is true for points on the
    # boundary of a domain, and false otherwise.

# The Dirichlet boundary values are defined using compiled expressions

# Define Dirichlet boundary (x = 0 or x = 1) expressions
c = Constant((0.0, 0.0, 0.0)) # More efficent to use Constant not Expression
    # Gamma_{D_0}: u = (0, 0, 0)
r = Expression(("scale*0.0",
                "scale*(y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta) - x[1])",
                "scale*(z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta) - x[2])"),
                scale = 0.5, y0 = 0.5, z0 = 0.5, theta = pi/3, degree=2)
    # Gamma_{D_1}: u=(0, (0.5+(y−0.5)cos(π/3)−(z−0.5)sin(π/3)−y)/2,
    #                    (0.5+(y−0.5)sin(π/3)+(z−0.5)cos(π/3)−z))/2)

# The boundary subdomains and the boundary condition expressions are collected
# together in two Dirichlet BC objects, one for each part of the boundary
bcl = DirichletBC(V, c, left)       # Gamma_{D_0}
bcr = DirichletBC(V, r, right)      # Gamma_{D_1}
    # The Dirichlet (essential) boundary conditions are constraints on the
    # function space V. The function space is therefore required as an argument
    # to DirichletBC
bcs = [bcl, bcr]                    # Combine all boundary conditions

# Trial and test functions, and the most recent approximate displacement, u are
# defined on the finite element space V.
# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
# Two constants are declared for the body force (B) and traction (T) terms
B  = Constant((0.0, -0.5, 0.0))  # Body force per unit volume
T  = Constant((0.1,  0.0, 0.0))  # Traction force on the boundary Gamma_N

    # In place of Constant it is also possible to use as_vector, e.g.
    # B = as_vector( [0.0,-0.5, 0.0] ). The advantage of Constant is that its
    # values can be changed without requiring re-generation and re-compilation
    # of C++ code. On the other hand, using as_vector can eliminate some
    # function calls during assembly.

# With the functions defined, the kinematic quantities involved in the model
# are defined using UFL syntax

# Kinematics
d = len(u)          # Length of displacement vector
I = Identity(d)     # Identity tensor
F = I + grad(u)     # Deformation gradient
C = F.T*F           # Right Cauchy-Green tensor C = F^T F

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Next, the material parameters are set and the strain energy density and the
# total potential energy are defined (using UFL syntax) 

# Elasticity parameters
E, nu = 10.0, 0.3                            # Young's Modulus, Poisson's ratio
mu = Constant(E/(2*(1 + nu)))                # Lame parameter
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu))) # Lame parameter

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

    # Just as for the body force and traction vectors, Constant has been used
    # for the model parameters 'mu' and 'lmbda' to avoid re-generation of C++
    # code when changing model parameters. Note that lambda is a reserved
    # keyword in Python, hence the misspelling 'lmbda'

# Directional derivatives are now computed of Pi and L (see equations
# 'first_variation' and 'second_variation')

# Compute first variation of Pi
F = derivative(Pi, u, v) # directional derivative about u in the direction of v
    # L(u;v) = dPI(u + epsilon nu)/ d(epsilon)
# Compute Jacobian of F
J = derivative(F, u, du)
    # a(u;du,v) = dL(u + epsilon du; v)/ d(epsilon)

# The complete variational problem can now be solved by a single call to solve
solve(F == 0, u, bcs, J=J) # L(u;v) == 0  a(u;du,v) = detF

# Finally, the solution 'u' is saved to a file named 'displacement.pvd' in VTK
# format, and the deformed mesh is plotted to the screen
file = File("displacement.pvd"); # Save solution in VTK format
file << u;

# Plot solution
plot(u)
plt.show()
