Overview

This repository contains MATLAB code for simulating Dictyostelium motion in clusters using 
the Cellular Potts Model(CPM) for cell motion and the Fitzhugh-Nagumo (FN) model for signal

Numerical Methods

-   Semi-implicit pseudo-spectrum method
-   Monte-Carlo


Components

1. Main.m - Main Simulation Script
2. Code_Excitationtime.m - Single cell FN simulation. Used for interpolating initial condition
3. energies_dicty.m - energy calculation for Cellular Potts Model
4. CPM_evlove.m - CPM evolution
5. Plot_speed - Code for plotting figures


To Successfully run the code:

Step 1 : Run Code_Excitationtime.m, which will create hilresult.mat recording the response of the FN model. 

Step 2 : Run Main.m, which is the main code and will create data for the cell center of mass
and the wave period.

Tune the number of cell in Main.m and get all data. 

The output includes
1. center of mass at different time points
2. spiral number 
3. wave period at different radius


After collecting all data, run "Plot_speed.m" to obtain figures for the angular velocity of waves and the angular velocity of cells.
Before plotting the figures, please create a folder with name "figure" to save all figures. Otherwise, you need to modified the save path in the code.
