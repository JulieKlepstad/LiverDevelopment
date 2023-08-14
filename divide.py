#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:03:20 2018

Author: Julie Klepstad
"""

import numpy as np
import random as rd


def divide(p, types, P_divide):
	N = len(p[:,0])	
	counter = 0
	attempts = 0
	while counter == 0:
		j = rd.choice(np.arange(0,N-1))
		rand = rd.uniform(0,1)
		if types[j] == 0 and rand < P_divide[0]:
			pj = p[j,:] +  np.random.uniform(-0.01,0.01, 3)
			p = np.vstack((p, pj))
			types = np.append(types, types[j])
			counter += 1
			
		if types[j] == 1 and rand > P_divide[0] and rand < (P_divide[0]+P_divide[1]):
			pj = p[j,:] +  np.random.uniform(-0.01,0.01, 3)
			p = np.vstack((p, pj))		
			types = np.append(types, types[j])
			counter += 1

		if types[j] == 2 and rand > (P_divide[0]+P_divide[1]):
			pj = p[j,:] +  np.random.uniform(-0.01,0.01, 3)
			p = np.vstack((p, pj))			
			types = np.append(types, types[j])
			counter += 1
			
		attempts += 1
		if attempts == N:
			break
		
	return p, types
