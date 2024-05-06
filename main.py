'''
                            MAIN.PY
    Main routine.       
    In this script you declare an evironment and initialize all_cells with elements of the Cell class.
    Then, the main loop is executed, which updates the enviroment every dt_enviroment, cell pheanomena
    every dt_phenotype and cell position and velocity every dt_mechanics.
    Every time the environment is updated, both the environment and the cells are showed in a plot. 
    These plots are saved in a folder named Figures.
    After, the total computation time is printed.
    Finally, results are saved: O2 field in the required points, number of cells and cell positions.
    If a video of the simulation is required, MakeVideo.py should be run.

'''

from fenics import *
from dolfin import *
import random
import matplotlib.pyplot as plt
import numpy as np
from cell import Cell
from constants import *
from variational_problem import Environment
import os
import time
from scipy.io import savemat

def plot_cells(E, all_cells, day, hour, minute):
    num_cells = len(all_cells)
    E.plot()
    plt.title(f"Day: {day}, Hour: {hour}, Min: {minute}, Cells: {num_cells}")
    positions = np.zeros((num_cells, 2))
    for n in range(num_cells):
        positions[n, :] = all_cells[n].position 
    plt.scatter(positions[:, 0], positions[:, 1], s = celular_radius, color='#FF0000',)

    plt.xlim([0, E.a])
    plt.ylim([0, E.b])    

    nombre_directorio = "Figures"
    if not os.path.exists(nombre_directorio):
        os.makedirs(nombre_directorio)
    os.chdir(nombre_directorio)
    cont = day*24*60+hour*60+minute
    plt.savefig(f"{cont}") 
    os.chdir("..")

    plt.show()

    return 


start_time = time.time()

# declare an environment.
E = Environment([width, height], nx, ny) 
E.setup()

# initialize the variable where the data about cells will be saved. 
# all_cells must be filled with elements of the Cell class.
all_cells = [] 
death_cells = 0
new_cells = 0

for i in range(num_initial_cells):
    cell = Cell(E)
    all_cells.append(cell)

# main loop
for minute in range(num_minutes):
    
    hour = (minute // 60) % 24
    minutes = minute % 60
    day = minute // (24 * 60)
    
    # O2 UPDATE (environment)
    if minute % dt_environment == 0:
        if not all_cells:
            print('ALL CELLS ARE DEAD')
            
        E.update(minute + 1, all_cells[0].uptake(E, all_cells))  # Update oxygen concentration
        plot_cells(E, all_cells, day, hour, minutes)
        
    # CELLS UPDATE
    
    if minute != 0 and minute % dt_mechanics == 0: # computation of cell velocity and position
        
        for cell in all_cells:
            cell.update_velocity(E, all_cells)
            cell.update_position(E, all_cells)           
    
    if minute != 0 and minute % dt_phenotype == 0: # phenomena related with phenotype: dead and proliferation
        
        for cell in all_cells:
            [death_cells, new_cells] = cell.update(E, all_cells, death_cells, new_cells)
   
end_time = time.time()
elapsed_time = end_time - start_time
elapsed_minutes = (elapsed_time // 60) % 60
elapsed_hours = (elapsed_time // 60) // 60
elapsed_seconds = elapsed_time % 60
print(f"The simulation took: {elapsed_hours:.2f} hours, {elapsed_minutes:.2f} minutes and {elapsed_seconds:.2f} seconds")
        

# Save results

# O2
filename = 'O2_fenics.mat'

yy = [25, 75, 125, 175, 225, 275, 325, 375] # points in y in which you want to save the oxygen
xx = [] # points in y in which you want to save the oxygen
for i_aux in range(47):
    xx.append(25*(i_aux+1))

irange = len(xx)
jrange = len(yy)
O2_fenics = np.zeros([jrange, irange])

for i in range(irange):
    
    for j in range(jrange):
        x = xx[i]
        y = yy[j]
        O2_fenics[jrange-1-j, i] = E.sol([x, y])
        
savemat(filename, {'O2': O2_fenics})

# Number of cells
filename = 'Cells.txt'
with open(filename, 'w') as file:
    ncells = len(all_cells) 
    file.write(str(ncells) + '\n') #Alive cells
    file.write(str(death_cells) + '\n') #Death cells
    file.write(str(new_cells) + '\n') #New cells

# Cell positions
filename = 'Positions.txt'

with open(filename, 'w') as file:
    for i in range(len(all_cells)):
        file.write(str(all_cells[i].position) + '\n')


