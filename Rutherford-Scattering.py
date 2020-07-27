#-----------------------------------
# Introduction
#-----------------------------------
#This program models Rutherford Scattering, i.e. particle scattering due to a central Kepler Force
#The procedure used is two-fold: first we model the motion of an incoming particle by solving Hamilton's Equations with an Euler-Method Algorithm
#Second, we determine the (approximate) polar scatering angle (between 0 and pi), and bin the result. 
#The result is the counts for each bin, which may be compared to the 1/sin**4(theta/2) relation.

#Written by Andrew Murphy, 2020

import numpy as np
import random
import matplotlib.pyplot as plt

#-----------------------------
# Functions used in Main Code
#-----------------------------

#Create a list of bins. Each bin is defined by the minimum value (bins[i][0]) and maximum value (bins[i][1])
def make_bins(bin_list):
    N = len(bin_list) - 1
    bins = []
    for i in range(N):
        bins.append([bin_list[i], bin_list[i + 1]])
    return bins

#Places a given value "val" into a bin
def bin_val(val, bins):
    for i in range(len(bins)):
        if bins[i][0] <= val < bins[i][1]:
            bin_count[i] += 1

# The following functions represent the six partial derivatives of the Hamiltonian for repulsive Kepler Motion: H = p**2/2m + k/r.  
def dH_dx(x,y,z):
    return -x*k/((x**2 + y**2 + z**2)**(1.5))

def dH_dy(x,y,z):
    return -y*k/((x**2 + y**2 + z**2)**(1.5))
    
def dH_dz(x,y,z):
    return -z*k/((x**2 + y**2 + z**2)**(1.5))
    
def dH_dpx(px, py, pz):
    return px/m

def dH_dpy(px, py, pz):
    return py/m

def dH_dpz(px, py, pz):
    return pz/m

#Here we define our scattering process: first we select a random start point within a beam-width of s_max, 
#then we approximate the motion of our particle by solving Hamilton's equations 
def scattering_process():
#Our first step is to randomly generate starting x and y coordinates within a circular cross section of radius s_max
    phi = random.uniform(-np.pi , np.pi)  #Random azimuthal angle
    r = random.uniform(0, s_max) #Random radius
    x = r*np.cos(phi) #x-coordinate of random point
    y = r*np.sin(phi) #y-coordinate of random point
    
    z = z_0
    px = 0
    py = 0
    pz = pz_0
    t = 0
    #Here we execute our Euler-Method Algorithm to approximately solve Hamilton's Equations in 3-Dimensions
    while t < t_max:
        diff_x = dH_dpx(px,py,pz)*dt
        diff_y = dH_dpy(px, py, pz)*dt
        diff_z = dH_dpz(px, py, pz)*dt
        diff_px = -dH_dx(x,y,z)*dt
        diff_py = -dH_dy(x,y,z)*dt
        diff_pz = -dH_dz(x,y,z)*dt
        
        x += diff_x
        y += diff_y
        z += diff_z
        px += diff_px
        py += diff_py
        pz += diff_pz
        
        t += dt
        
    return np.arccos(z/np.sqrt(x**2 + y**2 + z**2)) #Returns the polar angle of our particle at t_max, this is (approximately) our scattering angle
        
    

#--------------------------
# Seting up the Simulation
#--------------------------

particle_number = 1000

m = 1.0 # Particle Mass
k = 100.0 #Force-constant for Kepler-Force
z_0 = -100.0 #Initiapythonl Z-coordinate
pz_0 = 100.0 #Initial z-momentum
t_Steps = 1000 #Number of time steps
s_max = 1 #radius of particle beam

t_max = abs(4*z_0*m/pz_0) #This is our maximum time, which is sufficiently large, see explaination for why this is so
dt = t_max/t_Steps #This is our time-differential

#----------------------------
# Setting up the Bins
#----------------------------

number_of_bins = 100

bin_list = np.linspace(0, np.pi , (number_of_bins + 1))

bins = make_bins(bin_list)

bin_count = np.zeros(number_of_bins)

#-----------------------------
# Running the Simulation
#-----------------------------

for i in range(particle_number):
    bin_val(scattering_process(), bins)

#-----------------------------
# Plotting the Results
#-----------------------------

bin_coordinate = np.zeros(number_of_bins)

#Here we generate a list of the midpoint of each bin, the count of each bin is plotted at this point
for i in range(number_of_bins):
    bin_coordinate[i] = (bins[i][0]+ bins[i][1])/2
    
fig, ax = plt.subplots ()
ax.scatter(bin_coordinate, bin_count)
ax.set(xlabel='Scattering Angle (Rad)' , ylabel="Particle Count" , title='Scattering Results')
plt.show()
