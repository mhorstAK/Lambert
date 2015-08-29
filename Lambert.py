#!/usr/bin/env python
import math
import numpy as np
from math import atan2, pi, sqrt, sin, cos, cosh
from numpy import absolute, dot, array, sqrt, linalg, mod

r0 = [-96948447.3751163, 46106976.1901101, 6225649.83906096]
rf = [31161153.5750633, 143995536.213184, 3018.4216823707]
dt0 = 26179200
mu = 132712440018 #gravitational parameter of the sun
tol = .000001

#r0 = [31194146.7137977, 143988869.052709, 2352.10245497072]
#rf = [56033147.7402139, -781207927.044977, 1972452.68667642]
#dt0 = 103507200
#mu = 132712440018
#tol = .000001


mag_r0 = np.linalg.norm(r0)
mag_rf = np.linalg.norm(rf)
cdnu = dot(r0,rf)/(mag_r0*mag_rf)

nu1 = atan2(r0[1],r0[0])
nu2 = atan2(rf[1],rf[0])
dnu = nu2 - nu1
dnu = mod(dnu,2*pi)

if dnu < pi:
	DM = 1
else:
	DM = -1

A = DM*sqrt(mag_r0*mag_rf*(1+cdnu))

if dnu==0: 
	print('Delta Nu = 0!!!')

psi = 0
c2 = 0.5 
c3 = .16666666666667 
psi_up  = 4*(pi**2) 
psi_low = -4*pi
N = 0.8

eps = dt0 + 1000
while_count = 0

while eps > tol:
	y = mag_r0 + mag_rf + (A*(psi*c3 - 1))/sqrt(c2)

	#if np.any(A) > 0.0 and np.any(y) < 0.0:
	if A > 0 and y < 0:
		psi = N/c3 * (1 - sqrt(c2)/A * (mag_r0+mag_rf))
		
		if psi > tol:
			c2 = (1 - cos(sqrt(psi)))/psi
			c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi**3)
		elif psi < -tol:
			c2 = (1 - cosh(sqrt(-psi)))/psi
			c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)**3)
		else:
			c2 = 0.5
			c3 = .16666666666667

		y = mag_r0 + mag_rf + (A*(psi*c3 - 1))/sqrt(c2)

	chi = sqrt(y/c2)

	dt = ((chi**3)*c3 + A*sqrt(y))/sqrt(mu)

	if dt <= dt0:
		psi_low = psi
	else:
		psi_up = psi

	psi = 0.5*(psi_up + psi_low)

	if psi > tol:
		c2 = (1 - cos(sqrt(psi)))/psi
		c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi**3)
	elif psi < -tol:
		c2 = (1 - cosh(sqrt(-psi)))/psi
		c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)**3)
	else:
		c2 = 0.5
		c3 = .16666666666667

	eps = absolute(dt - dt0)
	while_count = while_count + 1

	if while_count > 1000:
		break

f = 1 - y/mag_r0
g_dot = 1 - y/mag_rf
g = A*sqrt(y/mu)

v0 = (np.array(rf) - np.multiply(f,r0))/g
vf = (np.multiply(g_dot,rf) - r0)/g

print(v0)
print(vf)
