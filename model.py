#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 7 14:56:15 2018

Author: Julie Klepstad
"""

# packages you need to download
import numpy as np # matrix operations
import copy as cp # make copies of arrays
import random as rd # makes random numbers
import os # to make folders
import matplotlib.pyplot as plt # to make plots
from mpl_toolkits.mplot3d import Axes3D # to make plots in 3D

# Import the other scripts. These have to be in the same folder as this model.py script
from neighbors import *
from move import *
from state import *
from divide import *
from startdata import *



#==============================================================================
# Settings for you to choose: Step 1/2
#==============================================================================
folder0 = 'model2/' # name of folder on the form 'name/'
N0 = 100 # initial number of cells
Nfolders = 10 # the number of simulations you want to run 



#==============================================================================
# Things you probably don't need to change:
#==============================================================================
if not os.path.exists(os.path.dirname(folder0)): # produces the folder if it does not exist
	os.makedirs(os.path.dirname(folder0), exist_ok=True)
	
# The final number of cells of each type in each simulation will be saved in a document cell_types.txt 	
with open(folder0 + 'cell_types.txt', 'w') as f: 
	f.write('Model number' +  '    ' +  'Number of timesteps' + '    ' +  'Number of progenitors' + '    ' + 'Number of H' + '    ' + 'Number of C' + '\n')

for iN in range(1,Nfolders+1):
	print('Starting running simulation number: ', iN)
	folder = folder0[0:-1] + '_' + str(iN) + '/'
	p = startdata(N0) # x,y,z, type (0: pregenitor, 1: H, 2: C), dt division
	p0 = cp.deepcopy(p)
	if not os.path.exists(os.path.dirname(folder)): # produces the folder if it does not exist
		os.makedirs(os.path.dirname(folder), exist_ok=True)	
	P_divide = np.zeros((3,1)) 
	dt = 0.1 # time resolution
	noise = 10**(-3) # order of magnitude of noise
	dt_save = 10	# how often we save the data
	T = 100000 # temporary number of timesteps
	counter_dt_save = 0
	l = np.ones((3, 3)) # strenghts of interactions: (affected cell, affecting cells)
	counter_dt_print = 0
	
	
	
	#==============================================================================
	# Initial parameters for you to choose: Step 2/2
	#==============================================================================
	N_phase2 = 200 # number of cells needed to change to phase 2
	Nmax = 3000 # maximum number of cells
	
	# probabilities of progenitor division outcomes
	P_divide[0] = 0.2 # progenitor becomes H and C
	P_divide[1] = 0.0 # progenitor becomes 2 C (BEC)
	P_divide[2] = 0.8 # progenitor becomes 2 H   
    
	# rates of dividing
	rate_progenitor = 1/100 
	rate_C = 1/100 
	rate_H = 1/100
    
	# set to True if the rates of dividing should change during the second phase
	change_division_times_during_simulation = False
	N_change_division = 2300
	rate_progenitor_2 = 1/100 
	rate_C_2 = 1/33 
	rate_H_2 = 1/100
	already_changed = False
	
	n = 10 # number of neighbors 
	n_nearest_neighbors = True # find the n closest cells
	
	# Wheter or not to plot a few images along the way to see the progression. Set to True for yes:
	plots = False 
	
	###############################################################################
	# Now you are done, the rest happens automaticly :)
	###############################################################################		
	
	
	
	
	
	#==============================================================================
	# Writes variables to a document that is saved in the file with the data
	#==============================================================================
	rates = np.array([rate_progenitor, rate_H, rate_C])
	rates_2 = np.array([rate_progenitor_2, rate_H_2, rate_C_2])
	with open(folder + 'variables' + '.txt', 'w') as f: 
		f.write('Number of cells at t0:' + '    ' + str(N0) + '\n')
		f.write('Number of cells at the beginning of phase 2:' + '    ' + str(N_phase2) + '\n')
		f.write('Number of cells at the end:' + '    ' + str(Nmax) + '\n')
		f.write('Time resolution, dt:' + '    ' + str(dt) + '\n')
		if n_nearest_neighbors == True:
			f.write('Neighbors are defined as: ' + str(n) + ' closest cells:' + '    ' + str(n_nearest_neighbors) + '\n')
		f.write('Order of magnitude of noise:' + '    ' + str(noise) + '\n')	
		f.write('Probability for progenitors to become H and C, 2 C, 2 H:' + '    ' + str(P_divide) + '\n')
		f.write('Rates of dividing for progenitors, C and H:' + '    ' + str(rates) + '\n')
		if change_division_times_during_simulation == True:
			f.write('Rates of dividing for progenitors, C and H after there are ' + str(N_change_division) + ' cells:' + '    ' + str(rates_2) + '\n')
		f.write('Strenght of interaction from p to p, H and C:' + '    ' + str(l[0,:]) + '\n')	
		f.write('Strenght of interaction from H to p, H and C:' + '    ' + str(l[1,:]) + '\n')	
		f.write('Strenght of interaction from C to p, H and C:' + '    ' + str(l[2,:]) + '\n')	
		
		
	#==============================================================================
	# Simulation
	#==============================================================================
	
	for it in range(0,T):
		# Find all cells' neighbors
		ID, types_n = neighbors_n(p, n, p[:,3]) # indexes of n closest neighbors, the types of all neighbors
		
		if n_nearest_neighbors == True:
			NB = cp.deepcopy(ID)
			NB_type = cp.deepcopy(types_n)
			av_neigh = n
		
		# Move according to the envrionment
		DR, V_average = move(p, l, p[:,3], NB, NB_type) # updates for positions, average potential per cell
		
		N = len(p[:,0])
		p[:,0:3] -= DR*dt + np.random.uniform(-noise,noise, (N,3)) # Euler
		
		if it == 0: # set times until first division
			p[:,4] = it*dt - rate_progenitor**(-1)*np.log(np.random.uniform(0,1,(N)))	
    		   
		if N < Nmax:
			if N < N_phase2: # progenitors divide and become progenitors
				for i in range(N):
					if p[i,4] >= it*dt and p[i,4] < (it+1)*dt:
						p[i,4] = 0 # reset clock
						pj = cp.deepcopy(p[i,:])
						pj[0:3] += np.random.uniform(-0.01,0.01, 3)
						p = np.vstack((p, pj))
		
			else:
				for i in range(N): # progenitors divide and become H and C, H and C divide and become them selfs
					if p[i,4] >= it*dt and p[i,4] < (it+1)*dt:
						p[i,4] = 0 # reset clock
						pj = cp.deepcopy(p[i,:])
						if p[i,3] == 0:
							rand = rd.uniform(0,1)
							if rand < P_divide[0]:
								p[i,3] = 1 # become H
								pj[3] = 2 # become C
							if rand > P_divide[0] and rand < P_divide[0]+P_divide[1]:
								p[i,3] = 2 # become C
								pj[3] = 2 # become C
							if rand > P_divide[0]+P_divide[1]:
								p[i,3] = 1 # become H
								pj[3] = 1 # become H							
								
						# if they are the same just move the copy and add it to the bulk		
						pj[0:3] += np.random.uniform(-0.01,0.01, 3)
						p = np.vstack((p, pj))
			
		# set a new division time for the cells that just divided
		N = len(p[:,0])
		for i in range(N): 
			if p[i,4] == 0:
				p[i,4] = it*dt - rates[int(p[i,3])]**(-1)*np.log(np.random.uniform(0,1))


		# change the rate of divisions
		if change_division_times_during_simulation == True and N >= N_change_division and already_changed == False:
			# update the rates of divisions for the cells if their cell type's rate of division has changed
			for i in range(N): 
				if p[i,3] == 0 and rates[0] != rates_2[0]: # progenitor cells
					p[i,4] = p[i,4] * rates[0]/rates_2[0]
				if p[i,3] == 1 and rates[1] != rates_2[1]: # H
					p[i,4] = p[i,4] * rates[1]/rates_2[1]
				if p[i,3] == 2 and rates[2] != rates_2[2]: # C
					print("time left until division: ", p[i,4]-it*dt, "new time until division :", (p[i,4]-it*dt) * rates[2]/rates_2[2])
					p[i,4] = it*dt + (p[i,4]-it*dt) * rates[2]/rates_2[2]
			already_changed = True
						
			rates = rates_2 # change the rates of divisions


		# save data
		counter_dt_save += 1
		if counter_dt_save == dt_save:
			t_string = "000000"
			t_string = t_string[:-len(str(it+1))]
			t_string = t_string + str(it+1)
			np.savetxt(folder + 'p_' + t_string + '.txt', p)
			counter_dt_save = 0
			
		counter_dt_print += 1
		if counter_dt_print == 1000:
			print(it, N, np.sum(p[:,3]==0),np.sum(p[:,3]==1),np.sum(p[:,3]==2), '\n')
			counter_dt_print = 0
		
		if N >= Nmax and counter_dt_save == 0:
			with open(folder0 + 'cell_types.txt', 'a') as f: 
				f.write(str(iN) +  '    ' +  str(it) + '    ' +  str(np.sum(p[:,3]==0)) + '    ' + str(np.sum(p[:,3]==1)) + '    ' + str(np.sum(p[:,3]==2)) + '\n')
			break


#==============================================================================
# Plot
#==============================================================================
	
if plots == True:
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1, projection='3d')
	index0 = np.argwhere(p[:,3]==0)
	index1 = np.argwhere(p[:,3]==1)
	index2 = np.argwhere(p[:,3]==2)
	ax.scatter(p[index0[:,0],0], p[index0[:,0],1], p[index0[:,0],2], s=100, color = 'blue')
	ax.scatter(p[index1[:,0],0], p[index1[:,0],1], p[index1[:,0],2], s=100, color = 'green')
	ax.scatter(p[index2[:,0],0], p[index2[:,0],1], p[index2[:,0],2], s=100, color = 'red')
	ax.legend(['Progenitor','H','C'])
