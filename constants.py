#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
                            CONSTANTS.PY
    Definition of the constants used in the different scripts of the code.

@author: mehran and lucia
"""
# Time steps and total time
dt_environment = 1 # Time step for environment
dt_mechanics = 1 # Time step for velocity and position
dt_phenotype = 6 # Time step for cells 

num_minutes = 1 * 24 * 60  # Number of minutes to run the simulation

# BC AND IC
u_0 = 7.0 # O2 initial condition (actually the conditions are defines in variational_problem) (mmHg)
u_l = 0.0 # Boundary conditions (left) (mmHg)
u_r = 7.0 # Boundary conditions (right) (mmHg)

# CONSTANS FOR VARIATIONAL_PROBLEM.PY
nx = 150 # Number of nodes of the mesh in the x direction
ny = 50 # Number of nodes of the mesh in the y direction  
width = 1200 # Width of the mesh (um)
height = 400 # Height of the mesh (um)
height_element = height/ny
width_element = width/nx
dz_physicell = 50 # Size of the element in the z direction (um)

nsteps = 10 # Steps used in the newton solver

# CONSTANS FOR CELL.PY
alpha_n = 60000 # mmHg.um³/cell.min
o2M = 2.5 # (mmHg) 
Diffusion_coefficient = 60000.0 

num_initial_cells = 250

celular_radius = 14.293/2 # (um)

birth_rate =  1/1440 # (1/min)
death_rate = 1/1440 # (1/min)
hipoxia_threshold =  7 # (mmHg)
anoxia_threshold = 1.6 # (mmHg)
sensitivity_to_anoxia = 0.1 # (mmHg)
equilibrium_spacing = 25.0 # (um) the new cells born appear more or less at this distance from the mothe cell
Kn = 45 #um2/(mmHg.min)
Dn = 3 #um2/min
persistence_time = 1 #min
