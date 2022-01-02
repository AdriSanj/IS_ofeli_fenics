import numpy as np
import dolfin
import matplotlib.pyplot as plt

# Setting matplotlib
plt.style.use('seaborn-poster')


k = 1.5 # wavenumber
distributed_load = dolfin.Expression("1.+0.*x[0]", degree=1) # linear polynomial expression
L = np.pi
n_cells = 10
n_nodes = n_cells + 1
mesh = dolfin.IntervalMesh(n_cells,0.,L)
plt.xlabel('$x$')
plt.title('finite element mesh')
h = dolfin.plot(mesh)
# Define the function space
V = dolfin.FunctionSpace(mesh,"CG", 1)

# Define the test and trial functions
v = dolfin.TestFunction(V)
u = dolfin.TrialFunction(V)
K_form = dolfin.inner(dolfin.grad(u), dolfin.grad(v)) * dolfin.dx
M_form = dolfin.inner(u, v) * dolfin.dx 

K = dolfin.PETScMatrix()
M = dolfin.PETScMatrix()

dolfin.assemble(K_form, tensor=K)
dolfin.assemble(M_form, tensor=M)

b = dolfin.Vector()

b_form = distributed_load * v * dolfin.dx
dolfin.assemble(b_form, tensor=b)

# Get discretization matrix
A = K - k**2*M

# Define and apply the boundary conditions 
bcs = dolfin.DirichletBC(V, dolfin.Constant(0), "on_boundary")
bcs.apply(A)
bcs.apply(b)

sol = dolfin.Function(V)
dolfin.solve(A, sol.vector(), b)
plt.figure(figsize=(10,8))
print(sol.vector().get_local())
dolfin.plot(sol)
x = mesh.coordinates()
plt.plot(x, (np.cos(k*(x-np.pi/2))/np.cos(k*np.pi/2) -1)/k**2, 'rs--')
plt.xlabel(r'$x$')
plt.ylabel(r'$u(x)$')
plt.legend(['Approx.','Exact'])
plt.show()
