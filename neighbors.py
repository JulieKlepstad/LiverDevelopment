#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 16:56:13 2018

@author: Julie Klepstad
"""

import numpy as np
import copy as cp

def neighbors_n(p, n, types):
	Farray = np.array
	Fwhere = np.where
	Fzeros = np.zeros
	Fsqrt = np.sqrt
	Fdelete = np.delete 
	
	# First, find the n number of closest cells:
	N = len(p[:,0])
	
	ID = Fzeros((N, n)) # neighbor array with indexes 
	types_n = Fzeros((N,n)) 

	r = Fsqrt( (Farray([p[:,0]]).T - p[:,0])**2 + (Farray([p[:,1]]).T - p[:,1])**2 + (Farray([p[:,2]]).T - p[:,2])**2 ) # N x N distance between all points

	# make the sorted matrix
	r_temp = cp.deepcopy(r)
	r_temp.sort(axis=1) # sorts all rows from small to large distances.
	r_temp = Fdelete(r_temp,0,1) # deletes first coloum = it self
	
	# Append all indexes of values below the inc'th value
	for ix in range(N):
		indexes = Farray(Fwhere(r[ix] < r_temp[ix,n])) # N x n+1 matrix with the closest cell distances. NOTICE that the cell itself is its own neighbor
		ID[ix] = indexes[indexes != ix] # fills in the indexes without including it self
		for iy in range(n):
			if types[int(ID[ix, iy])] == 0:
				types_n[ix, iy] = 0
			elif types[int(ID[ix, iy])] == 1:
				types_n[ix, iy] = 1
			else:
				types_n[ix, iy] = 2

	ID = ID.astype(int) # floats to ints

	return ID, types_n
	
