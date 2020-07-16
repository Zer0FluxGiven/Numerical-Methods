#-----------------------------------------
#   Introduction
#-----------------------------------------
#This program attempts to model magnetic ordering in materials (e.g. Ferromagnetism and Anti-Ferromagnetism) with the 2D Ising Model
#The model consists of an NxN square lattice, with each point having a spin attributed to it (spin-up = +1, spin-down = -1)
#The spins of each lattice point interact with their nearest neighbor, as well as with an external applied magnetic field (H-Field)

#Written by Andrew Murphy, 2020

import numpy as np
import random
import matplotlib.pyplot as plt

N = 10 #Here we set our ensemble size, which will be an NxN square lattice

#Here we define our interaction parameter J, 
#when J = 1, we expect the model to exhibit Ferromagnetic behavior
#when J = -1, we expect the model to exhibit Anti-Ferromagnetic behavior
J = 1

#----------------------------------
#   Functions used in the Model
#----------------------------------

#Here we generate our NxN ensemble lattice. Initially, the spins are randomly oriented
def make_ensemble():
    ensemble = [[0 for y in range(N)] for x in range(N)]
    i = 0
    while i < N :
        j = 0
        while j < N :
            ensemble[i][j] = random.choice([-1, 1])
            j = j+1
        i = i+1
    return ensemble
    
# Here we define our interaction Hamiltonian for a lattice point (i,j)
def H_int(i,j): 
	return -J*ensemble[i][j]*(ensemble[(i+1)%N][j] + ensemble[(i-1)%N][j] + ensemble[i][(j+1)%N] + ensemble[i][(j-1)%N])

#Here we define our Metropolis Algorithm, 
#this algorithm randomly selects a lattice point and determines whether or not the spin is flipped
def mc_metropolis():
    i = random.randint(0,N-1) #x-coordinate of random lattice point
    j = random.randint(0,N-1) #y-coordinate of random lattice point
    
    microstate = ensemble[i][j]
    
    E_0 = (-H_field * microstate) + H_int(i,j) # Energy of unflipped state
    E_f = (H_field * microstate) - H_int(i,j) # Energy of flipped state
    dE = E_f - E_0
    
    #This is the crux of the algorithm: 
    #if the energy is lowered by flipping, the spin is flipped
    #otherwise, the spin flips with probability exp(-dE/T)
    if dE < 0 :
        ensemble[i][j] = -microstate
    elif random.random() < np.exp(-dE/T) :
        ensemble[i][j] = -microstate
    return ensemble
 
#Here we define a function to calculate the energy of the ensemble by summing over each point
def calculate_energy():
    energy = 0
    i = 0
    while i < N :
        j = 0
        while j < N :
            energy = energy + (-H_field * ensemble[i][j]) + ((0.25)*H_int(i,j))
            j = j+1
        i = i +1 
    return energy
    
# Here we define a function to calculate the magnetization of the ensemble    
def calculate_magnetization():
    magnetization = 0
    i = 0
    while i < N :
        magnetization = magnetization + (sum(ensemble[i])/(N**2))
        i = i + 1
    return magnetization
    
#-------------------------------
#   Main Part of Code
#-------------------------------

#Here we set up the parameters of the model:
# Temperature (T) and External Magnetic Field (H-Field) 
# as well as the sweeping parameters for Temp-Sweeps and Field-sweeps

T = 0.1 #Be sure, T =/= 0 !
T_step = 0.01
T_max = 5.1

H_field = 0
H_field_step = 0.1
H_field_max = 5

#These values determine whether the model performs a Temperature-Sweep or a Field-Sweep
Temp_Sweep = 1 # =1 for Temperature Sweep
Field_Sweep = 0 #=1 for Field Sweep

#These values determine the number of Monte-Carlo steps taken for each step of the sweep (temp or field
eqSteps = 10000
measureSteps = 10000

T_vals = [] #Temperature 
M_vals = [] #Magnetization
C_vals = [] #Heat Capacity
H_vals = [] #H-Field
E_vals = [] #Internal Energy

ensemble = make_ensemble()

if Temp_Sweep == 1 :
    while T < T_max:
        E1 = E2 = M = 0
        eq = 0
        while eq < eqSteps : #Equillibriate the ensemble
            ensemble = mc_metropolis()
            eq = eq + 1
        measure = 0
        while measure < measureSteps: #Conduct measurements of Heat Capacity and Magnetization
            ensemble = mc_metropolis()
            energy = calculate_energy()
            mag = calculate_magnetization()
            E1 = E1 + energy
            E2 = E2 + energy**2
            M = M + mag
            measure = measure + 1
        avg_energy = E1/measureSteps
        avg_mag = M/measureSteps
        capacity = ((E2/measureSteps) - ((E1**2)/(measureSteps**2)))/(T**2)
        M_vals.append(avg_mag)
        C_vals.append(capacity)
        T_vals.append(T)
        
        T = T + T_step
    
if Field_Sweep == 1 :
    while H_field < H_field_max:
	E1 = E2 = M = 0
        eq = 0
        while eq < eqSteps : #Equillibriate the ensemble
            ensemble = mc_metropolis()
            eq = eq + 1
        measure = 0
        while measure < measureSteps: #Conduct measurements of Heat Capacity and Magnetization
            ensemble = mc_metropolis()
            energy = calculate_energy()
            mag = calculate_magnetization()
            E1 = E1 + energy
            E2 = E2 + energy**2
            M = M + mag
            measure = measure + 1
        avg_energy = E1/measureSteps
        avg_mag = M/measureSteps
        capacity = ((E2/measureSteps) - ((E1**2)/(measureSteps**2)))/(T**2)
	#Here we calculate the heat capacity, using the equation C = (stdev(E)/T)**2
        M_vals.append(avg_mag)
        C_vals.append(capacity)
        H_vals.append(H_field)
        
        H_field = H_field + H_field_step
   
#Here we plot our values. In my example, I simply plot the temperature and heat capacity.
plt.scatter(T_vals, C_vals)
plt.show()
