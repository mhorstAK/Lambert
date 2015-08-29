#!/usr/bin/env python
import math
import numpy as np
from math import atan2, pi, sqrt, sin, cos, cosh
from numpy import absolute, dot, array, sqrt, linalg, mod

class Lambert:
	"""docstring for Lambert"""
	#Listed are variables which are consistant with in both of the two cases
	mu = 132712440018 #gravitational parameter of the sun
	tol = .000001 #tolerance thingy (I think)

	psi = 0
	c2 = 0.5 
	c3 = .16666666666667 
	psi_up  = 4*(pi**2) 
	psi_low = -4*pi
	N = 0.8

	def __init__(self, r0, rf, dt0):
	#	super(Lambert, self).__init__()
		"""Creates Vaiables Which Change Per Case"""
		self.r0 = r0
		self.rf = rf
		self.dt0 = dt0

	def mag_f(self):
		"""Normalizes the Vectors"""
		mag_r0 = np.linalg.norm(self.r0)
		mag_rf = np.linalg.norm(self.rf)
		print(mag_r0)
		print(mag_rf)

	#def cdnu_f(self, mag_f):
	#	"""Lets be real, nobody gets this."""
	#	cdnu = dot(self.r0,self.rf)/(self.mag_r0*self.mag_rf)

#r0 = [-96948447.3751163, 46106976.1901101, 6225649.83906096]
#rf = [31161153.5750633, 143995536.213184, 3018.4216823707]
#dt0 = 26179200

#lam = Lambert()
#lam.mag_f([-96948447.3751163, 46106976.1901101, 6225649.83906096], [31161153.5750633, 143995536.213184, 3018.4216823707], 26179200)


