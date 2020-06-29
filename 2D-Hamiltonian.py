#This program solves Hamilton's Equations in 2 Degrees of Freedom, and plots the phase space trajectory in two 2-Dimensional plots.
#The program uses a Euler-Method Algorithm to generate the position and momentum coordinates of a particle's trajectory in phase space

# Written by Andrew Murphy, 2020

import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def dH_dx1(x1,x2, p1, p2) :                #Here we define the partial derivatives of the Hamiltonian with repsect to the two position 			
    return 0                               #and two momentum variables. In this example, we use the Hamiltonian                                 
                                           #H = (p1**2 + p2**2)/2m + mg*x2 for projectile motion
def dH_dx2(x1,x2, p1, p2) :            
    m = 0.5
    g = 9.81
    return  m*g
	
def dH_dp1(x1,x2, p1, p2) :            
    m = 1                 			         			
    return p1/m
	
def dH_dp2(x1,x2, p1, p2) :            
    m = 1                			         			
    return p2/m

    
	
x1_0 = 0
x2_0 = 0             # These are our inital parameters, initial position-coordinates (x1_0,  x2_0),
p1_0 = 15
p2_0 = 15         # initial momenmtum-coordinates (p1_0, p2_0), maximum time (tmax), and number of iterations (n_max) 
tmax  = 5          # Note: initial time is assumed to be 0
n_max = 10000

x1val = [x1_0]
x2val = [x2_0]
p1val = [p1_0]
p2val = [p2_0]

x1 = x1_0
x2 = x2_0             #Here we set up our variables that we will use in our Euler-Method Algorithm
p1 = p1_0
p2 = p2_0
n = 0
dt = tmax / n_max 

while n < n_max :
    x1_new = x1 + dH_dp1(x1, x2, p1, p2) * dt   #Here we execute our Euler-Metohd Algorithm, using Hamilton's equations:
    x2_new = x2 + dH_dp2(x1, x2, p1, p2) * dt   # dxi/dt = dH/dpi ; dp/dt = - dH/dxi 
    p1_new = p1 - dH_dx1(x1, x2, p1, p2) * dt
    p2_new = p2 - dH_dx2(x1, x2, p1, p2) * dt
    
    x1 = x1_new
    x2 = x2_new
    p1 = p1_new
    p2 = p2_new
    
    x1val.append(x1)
    x2val.append(x2)
    p1val.append(p1)
    p2val.append(p2)
    n = n + 1	

fig, axs = plt.subplots(2)          #Here we plot the trajectories through phase-space in two plots, stacked vertically
axs[0].plot(x1val, x2val)           # axs[0] (top plot) is the position trajectory (x1 on the x-axis, x2 on the y-axis)
axs[0].set_title("Position")        # axs[1] (bottom plot) is the momentum trajectory (p1 on the x-axis, p2 on the y-axis)
axs[1].plot (p1val, p2val)
axs[1].set_title("Momentum")
plt.show()            
