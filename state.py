# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:23:37 2018

Author: Julie Klepstad
"""


import numpy as np
import random as rd

def state(p, types, P_H, random_switch, neighbor_switch, NB, NB_type):
	counter = 0
	if random_switch == True:
		if np.any(types == 0): # change state of 1 progenitor cell, if there are any left
			while counter == 0:
				NB0 = np.argwhere(types == 0)
				j = rd.choice(NB0[:,0])
				rand = rd.uniform(0,1)
				if types[j] == 0:
					if rand < P_H:
						types[j] = 1
					else:
						types[j] = 2
					counter += 1

	if neighbor_switch == True:
		if np.any(types == 0): # change state of 1 progenitor cell, if there are any left
			while counter == 0:
				NB0 = np.argwhere(types == 0)
				j = rd.choice(NB0[:,0])
				if types[j] == 0:
					NBj = NB[j][NB[j] >= 0] # the indexes of the neighbors of deep cell i
					N_NB = len(NBj)
					NBj_type = NB_type[j][0:N_NB] # the type of the neighbors of deep cell i
					NB1 = np.sum(NBj_type == 1)
					NB2 = np.sum(NBj_type == 2)
					N_non_progenitor = NB1 + NB2
					if N_non_progenitor > 0:
						if NB1/N_non_progenitor < P_H:
							types[j] = 1
						else:
							types[j] = 2
						counter += 1
					else:
						rand = rd.uniform(0,1)
						if rand < P_H:
							types[j] = 1
						else:
							types[j] = 2
						counter += 1						
	return types

