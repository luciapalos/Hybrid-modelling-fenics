'''
                            CELL.PY
    Script that initializes cells with its main features.
    This script also computes everything related to cells:
        - Phenomena such as death and proliferation
        - Velocity and position update
        - Cell uptake (this is then used in variational_problem.py)

'''

import random
from math import sqrt
from constants import *
from fenics import *
from variational_problem import Environment
import numpy as np
import copy

from scipy.spatial import KDTree
import matplotlib.pyplot as plt

class Cell:
    def __init__(self, E):
 
        pos_x = 0 + np.random.uniform()*400 
        pos_y = 0 + np.random.uniform()*400
        self.position = np.array([pos_x, pos_y])
            
        self.velocity = np.array([0.0, 0.0])
        self.equilibrium_spacing = equilibrium_spacing
        self.birth_rate = birth_rate
        self.death_rate = death_rate
        self.anoxia_threshold = anoxia_threshold
        self.hipoxia_threshold = hipoxia_threshold
        self.sensitivity_to_anoxia = sensitivity_to_anoxia
        self.dt = dt_phenotype
        self.dt_v = dt_mechanics
        self.persistence_time = persistence_time
        
        self.Kn = Kn
        self.Dn = Dn


    def update_velocity(self, E, all_cells):
        
        o2 = E.sol(self.position)
        
        # F_loc
        a = random.random()
        
        if (a <= self.dt_v/self.persistence_time):
                    
            #RANDOM MIGRATION
            b = 0.5 # parameter that graduates randomness in migration
            vrandom_aux = np.random.uniform(-1,1,2) # new direction (random):
            vrandom = vrandom_aux/np.linalg.norm(vrandom_aux)
            
            #CHEMOTAXIS
            mod_o2 = np.linalg.norm(E.gradu(self.position))
            d = E.gradu(self.position)/mod_o2
    
            #ULOC 
            uloc_aux = b*d+vrandom*(1-b)
            uloc = uloc_aux/np.linalg.norm(uloc_aux)
            
            if o2 < self.hipoxia_threshold:
                Xn = (1-o2/self.hipoxia_threshold)
            elif o2 >= self.hipoxia_threshold: 
                Xn = 0.0
                               
            sloc = 0.5
            v_loc = sloc*uloc
 
            self.velocity = v_loc 
            
            
        # BOUNDARY CONDITIONS: We have to take into account that the cells cannot go out of the domain
        
        dist_xr = E.a - self.position[0] # distance to the right boundary
        dist_xl = self.position[0] # distance to the left boundary
        dist_yu = E.b - self.position[1] #distance to the upper boundary
        dist_ylow = self.position[1] #distance to the low boundary 
        
        if np.linalg.norm(self.velocity) > (celular_radius/self.dt_v): # if, at this velocity, the cell covers a distance which is greater than the celular radius
            
            # Check if the distance to the border is greater than the distance covered in a dt_v
            if ((dist_xr < self.dt_v*self.velocity[0]) and self.velocity[0]>0) \
               or ((dist_xl < abs(self.dt_v*self.velocity[0])) and self.velocity[0]<0):
                self.velocity[0] = -self.velocity[0]
            
            if ((dist_yu < self.dt_v*self.velocity[1]) and self.velocity[1]>0)\
                or (((dist_ylow < abs(self.dt_v*self.velocity[1]))) and self.velocity[1]<0):
                self.velocity[1] = -self.velocity[1]

        else:
            # Check if the distance to the border is greater than the celular radius
            if ((dist_xr < celular_radius) and self.velocity[0]>0) or ((dist_xl < celular_radius) and self.velocity[0]<0):
                self.velocity[0] = -self.velocity[0]
                
            if ((dist_yu < celular_radius) and self.velocity[1]>0) or ( (dist_ylow < celular_radius) and self.velocity[1]<0):
                self.velocity[1] = -self.velocity[1]
                         

    def division(self, all_cells, E): # function that regulates proliferation

        c = copy.deepcopy(self)
        all_cells.append(c)
        r = self.equilibrium_spacing
        angle = np.random.uniform()
        c.position[0] += r * np.cos(angle)
        c.position[1] += r * np.sin(angle)
        
        # loop to avoid the new cell appearing outside the domain
        if c.position[0]<0:
            c.position[0] = c.position[0] - c.position[0] + r * np.cos(angle)
        elif c.position[0]>E.a :
            c.position[0] = c.position[0] - (c.position[0] - E.a) - r * np.cos(angle)
            
        if c.position[1]<0:
            c.position[1] = c.position[1] - c.position[1] + r * np.sin(angle)
        elif c.position[1]>E.b :
            c.position[1] = c.position[1] - (c.position[1] - E.b) - r * np.sin(angle)
        
    
    def update(self, E, all_cells, death_cells, new_cells):
        
        o2 = E.sol(self.position)
        
        # PROLIFERATION
        hipoxia = o2<self.hipoxia_threshold
        if hipoxia:
            birth_rate = self.birth_rate *o2/self.hipoxia_threshold
        elif not hipoxia:
            birth_rate = self.birth_rate
           
        if birth_rate < 0:
            birth_rate = 0.0
                
        prob_birth = birth_rate * self.dt
        birth_decision = np.random.uniform(size=1)
        
        if (birth_decision <= prob_birth):
            self.division(all_cells, E)
            new_cells = new_cells + 1
            
        # DEATH
        death_rate = self.death_rate*0.5*(1-np.tanh((o2-self.anoxia_threshold)/self.sensitivity_to_anoxia))
        tangent = 1-np.tanh((o2-self.anoxia_threshold))
        
        if death_rate < 0.0:
            death_rate = 0.0
        
        prob_death = death_rate * self.dt
    
        death_decision = np.random.uniform(size=1)[0]

        if (death_decision <= prob_death):
            all_cells.remove(self) # if the cell dies, it is removed from all_cells
            death_cells = death_cells + 1
   
        return death_cells, new_cells
  
           
    def update_position(self, E, all_cells): 
        self.position = self.position + self.dt_v  * self.velocity

    # Function that regulates O2 uptake. It takes the O2 function from variational_problem.py. With this and the number of cells, it gives back the uptake to variational_problem.py.    
    def uptake(self, E, all_cells): 
        
        # Initialization of the arrays which will store the uptake and number of cells associated to each vertex. 
        uptake_values = np.zeros(E.mesh.num_vertices())
        num_cells = np.zeros(E.mesh.num_vertices())
        
        # Loop that for each vertex computes the cells associated it (the cells that are nearer from it that from other vertex).
        # With the number of cells, the uptake is computed. 
        # The values are stored in a way that FEniCS can understand. For more infomration, see the difference between degrees of freedom and vertices in FEniCS.
        for vertex in vertices(E.mesh):
            node_x = vertex.point().x()
            node_y = vertex.point().y()  
            
            for cell in all_cells:
                cell_x, cell_y = cell.position
                distance_x = abs(cell_x - node_x)
                distance_y = abs(cell_y - node_y)
                           
                if (distance_x < width_element/2) and (distance_y < height_element/2):
                    num_cells[E.v2d[vertex.index()]] += 1
                    
            cell_count = num_cells[E.v2d[vertex.index()]]
            
            concentration = E.C[E.v2d[vertex.index()]]
            uptake_coefficient = alpha_n * concentration/ (o2M + concentration) * cell_count/(width_element*height_element*dz_physicell)
            uptake_values[E.v2d[vertex.index()]] = uptake_coefficient
    
        return uptake_values
    