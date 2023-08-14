#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 20:56:54 2018

author: Julie Klepstad
"""

import numpy as np

def move(p, l, types, NB, NB_type):

	Fexp = np.exp
	Fsqrt = np.sqrt
	Fzeros = np.zeros 
	Farray = np.array
	Fmean = np.mean
	
	beta = 5 # range of attraction
	N = len(p[:,0])	
	DR = Fzeros((N, 3)) # for change in x and y and z
	V = Fzeros((N)) # potential
	
	for i in range(N):
		Pi = p[i,0:3]
		NBi = NB[i][NB[i] >= 0] # the indexes of the neighbors of i

		type_j = NB_type[i][0:len(NBi)] 

		xyzj = p[NBi,0:3] # positions of all neighbors
		xyz = Farray(xyzj - Pi[0:3])
		d = Fsqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2)
		xyzrel = xyz/d.reshape((len(d),1)) # relative distances

		# calculate derivative of V (sum of V from all neighbors):
		derivative_V = Fzeros((3))

		for k in range(len(NBi)):
			der_factor1 = Fexp(-d[k]/beta)

			l_np = l[int(types[i]), int(type_j[k])]
			
			derivative_V += der_factor1 * ((Fexp(-(beta-1)*d[k]/beta) - l_np/beta)*xyzrel[k])
			V[i] += Fexp(-d[k]) - l_np*der_factor1
	
		DR[i] = derivative_V		
		V_average = Fmean(V)

	return DR, V_average

	
