# This program solves Hamilton's Equations in 1D using an Euler-Method Algorithm
# The results of the algorithm ( xval and pval ) are phase space coordinates of position and momentum
# We plot the solution with position on the x-axis and momentum on the y-axis

# Written by Andrew Murphy, 2020

import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np

def dH_dx(x,p) :            # Here we are defings the partial derivatives of the Hamiltonian
    m = 0.5                 # In this case, we have a Harmnic Oscillator, H = p**2 / 2m + mw**2/2 * x**2
    w = 2 * np.pi           # m and w are mass and angfrequnecy 
    return m * (w**2) * x
	
def dH_dp(x,p) :
    m = 0.5
    return p/m
    
	
x_0 = 0             # These are our inital parameters, initial x-coordinate (x_0),
p_0 = 1.5           # initial p-coordinate (p_0), maximum time (tmax), and number of iterations (n_max) 
tmax  = 10          # Note: initial time is assumed to be 0
n_max = 10000

xval = [x_0]
pval = [p_0]

x = x_0             #Here we set up our variables that we will use in our Euler-Method Algorithm
p = p_0
n = 0
dt = tmax / n_max 

while n < n_max :
	x_new = x + dH_dp(x,p) * dt     # Here we execute our Euler-Method Algorithm, using Hamilton's Equations
	p_new = p - dH_dx(x, p) * dt    # dp/dt = -dH/dx and dx/dt = dH/dp
	
	x = x_new
	p = p_new
	xval.append(x)
	pval.append(p)
	n = n + 1	

plt.plot(xval,pval)     # Here we plot the solution, position is represented on the x axis, 
plt.show()              # And momentum is on the y-axis 
