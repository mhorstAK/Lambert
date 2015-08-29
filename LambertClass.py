#!/usr/bin/env python
import math
import numpy as np
from math import atan2, pi, sqrt, sin, cos
from numpy import absolute, dot, array, sqrt, mod

class Lambert:
	def __init__(self, r0, rf, dt0, tol):
		self.r0 = r0
		self.rf = rf
		self.dt0 = dt0
		self.tol = tol

	def Initialize_Solver(self, r0, rf):
		""" Normalizing vectors, finding the difference between vectors
		...
		"""
		mag_r0 = sqrt(dot(r0[0], r0[1]))
		mag_rf = sqrt(dot(rf[0], rf[1]))
		cosDV = dot(r0,rf)/(mag_r0*mag_rf)

		v1 = atan2(r0[0],r0[1]) #rad
		v2 = atan2(rf[0],rf[1]) #rad
		dv = v2 - v1
		dv = mod(dv, 2*pi)

		return mag_r0, mag_rf, cosDV, dv
