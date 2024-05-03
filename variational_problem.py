'''
                            VARIATIONAL_PROBLEM.PY
    In this script, the Environment class is initialized.
    A setup is done, in which the main features of the variational problem are defined: mesh, function space, bc, ic...
    Then, the function update, updates de O2 field.
    Finally with the function plot you can see the O2 field in a plot.

'''

import numpy as np
import matplotlib.pyplot as plt
from fenics import *
import copy
from matplotlib.pyplot import colorbar
from constants import *

# Declare environment class
class Environment:
    def __init__(self, shape=[width, height], m=nx, n=ny):
        self.a = shape[0]
        self.b = shape[1]
        self.m = m
        self.n = n
        self.t = 0.0
        self.dt = dt_environment/nsteps
        self.C = np.zeros((self.m+1) * (self.n+1))
        self.sol = []
        self.gradu = []
        self.D = Diffusion_coefficient

        return

    def setup(self):
        
        # Define the mesh
        mesh = RectangleMesh(Point(0, 0), Point(self.a, self.b), self.m, self.n)
        
        # Define the scalar function space (used for the solution) and the vector function space (used for the gradient of the solution)
        self.W_D1 = FunctionSpace(mesh, 'P', 1)
        self.Vgrad = VectorFunctionSpace(mesh, 'P', 1)
        
        # Define the TestFunction and the Function of the variational problem
        self.r = TestFunction(self.W_D1)
        self.o2 = Function(self.W_D1)
        self.o2_old = Function(self.W_D1)
        
        self.mesh = mesh
        
        # Functions used to map vertex and nodes in FEniCS
        self.v2d = vertex_to_dof_map(self.W_D1)
        self.d2v = dof_to_vertex_map(self.W_D1)
        
        # Define the boundary conditions
        width = self.a
        def boundary_left(x, on_boundary):
            return near(x[0], 0)                
        def boundary_right(x, on_boundary):
            return (near(x[0], width))
    
        u_left = Constant(u_l)
        u_right = Constant(u_r)
        bc_left = DirichletBC(self.W_D1, u_left, boundary_left)
        bc_right = DirichletBC(self.W_D1, u_right, boundary_right)
        self.bcs_o2 = [bc_left, bc_right]
                
        # Define the initial conditions
        u_aux = Expression('7*x[0]/L',degree=1, L = self.a)
        u_0 = interpolate(u_aux,self.W_D1)
        
        self.o2.assign(u_0)
        self.C = self.o2.vector()[:]
        self.o2_old.assign(u_0)
        
        return

    def update(self, t_o2, uptake):

        # Define a Function to store the uptake values
        uptake_function = Function(self.W_D1)
        uptake_function.vector()[:] = uptake
        
        # Define the fuction you want to solve
        self.F_D1 = (
            inner(self.o2 - self.o2_old, self.r) * dx
            + self.dt * self.D * inner(grad(self.o2), grad(self.r)) * dx
            + self.dt * uptake_function * self.r * dx
        )

        # Solver
        for time in np.arange(0, dt_environment, self.dt):
            set_log_level(30) # to avoid showing all the iterations
            solver_parameters = {
                "newton_solver": {
                    "linear_solver": "gmres",
                    "preconditioner": "ilu",
                    "maximum_iterations": 1000,
                    "relaxation_parameter": 0.7,
                    "relative_tolerance": 1e-7,
                    "absolute_tolerance": 1e-8
                }
            }
            
            solve(self.F_D1 == 0, self.o2, self.bcs_o2, solver_parameters=solver_parameters)  
            
            self.C = self.o2.vector()[:] # store the solution in a numpy array
            self.C[self.C<0.0] = 0.0 # avoid negative values
            
            self.sol = self.o2
            
            self.t += self.dt
            self.o2_old.assign(self.o2)

            self.gradu = project(grad(self.sol),self.Vgrad) # compute the gradient

    def plot(self): # function to plot the environment
        plt.figure(2)
        p = plot(self.o2, show=False)
        plt.colorbar(p)