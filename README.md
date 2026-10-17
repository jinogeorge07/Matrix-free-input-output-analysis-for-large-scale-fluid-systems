# gmres matrix free resolvent
Gmres matrix-free adjoint based method for input-output analysis

The main code is in the folder "Code for wavy wall and channel flow"

The code "main_2D_resolvent_precondition_1D_resolvent_wavywall_final.m" will run the wavy wall case
The code "main_2D_xy_resolvent_precondition_1D_resolvent_final.m" will run the channel flow case

This folder contains required subroutines, including Chebyshev and Fourier differentiation matrices. This folder also contains the mean flow required as input for wavy wall case. 

Each folder with the name "figure X" will contain source data (in .mat format), and post-processing code that convert this to .csv file and generate the figure. 