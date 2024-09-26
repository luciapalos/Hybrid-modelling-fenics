## Hybrid-modelling-fenics
Hybrid code: cancer cells are modelled with a discret model and oxygen is simulated with a continuum model.

# main.py
Script that has to be run in order to start the whole code.
In this script, you declare an environment (imported from variational_proble.py), an object that in this case represents the oxygen concentration in the domain. Then, the list all_cells is initialized, in which objects of Cell (imported from cell.py) type are saved. 
After, the main loop starts. In this loop, the oxygen is updated (every dt_environment), calling the function update from the object Environment. In order to do this, the oxygen uptake has to be computed, which is done by calling all_cells[0].uptake(E,all_cells). Then, every dt_mechanics cell velocity and position is updated and every dt_phenotype cell proliferation and death is updated.
The results are plotted every minute.
When the main loop ends, theres is a block of code to save the results:
-O2 is saved in a .mat file (the points where O2 is saved can be chosen modifying the arrays yy and xx).
-The number of alive, death and total cells is saved in a .txt file.
-The poistion of cells is saved in a .txt file.

# variational_problem.py
In this script, the Environment class is defined. It will be used to solve O2 with a FEM.
When the object is first called, the attributes of the Environment are initialized: size(a,b), number of elements in the mesh (m,n) and several lists to save results.
A setup is done, in which the main features of the variational problem are defined: mesh, function space, bc, ic...
The functions update is used to update the O2 concentrations every time that is called from main.py. A solver is used for this.
With the funcion plot you can serr the O2 field in a plot.

# cell.py
In this script, the Cell object is defined. 
When first called, the main variables of the object are initialized: positions, velocity, death and birth rates, etc.
Then, when the function update is called, the state of the cell (death or divide) is updated. 
When update_velocity or update_position are called, the velocity and position are updated respectively.
The function uptake is called from main.py and computes the oxygen uptake in each node. This is an input for variational_problem.py, because the oxygen consumption is needed to solve the oxygen evolution.

# constants.py
In this script the constants used in the rest of the code are defined.

# MakeVideo.py
This is an independent script (which means that it will not run when you run main.py) that you should run when main.py finish in order to get a video with the different plots of the simulation (this plots will save in a folder name Figure, this can be seen in main.py).
