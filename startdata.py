#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 12:45:00 2018

Author: Julie Klepstad
"""

import numpy as np

def startdata(N):
	p = np.zeros((N,5)) # x,y,z, type, dt division
	
	N = len(p[:,0])
	phi = np.random.uniform(0, 2*np.pi, N)
	costheta = np.random.uniform(-1, 1, N)
	u = np.random.uniform(0, 1, N)
	
	theta = np.arccos(costheta)
	r = 3 * np.power(u, 1/3)
	
	p[:,0]  = r * np.sin(theta) * np.cos(phi)
	p[:,1]  = r * np.sin(theta) * np.sin(phi)
	p[:,2]  = r * np.cos(theta)
	
	return p