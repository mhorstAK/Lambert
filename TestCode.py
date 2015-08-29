#!/usr/bin/env python
import math
import numpy as np
from math import atan2, pi, sqrt, sin, cos, exp, sinh
from numpy import absolute, dot, array, mod

r0 = [50,70]
rf = [40,80]
tol = 1e-6
dt0 = 153
mu = 2

mag_r0 = sqrt(dot(r0[1], r0[0]))
mag_rf = sqrt(dot(rf[1], rf[0]))
cosDV = dot(r0,rf)/(mag_r0*mag_rf)
print(cosDV)

v1 = atan2(r0[1],r0[0]) #rad
v2 = atan2(rf[1],rf[0]) #rad
dv = v2 - v1
dv = mod(dv, 2*pi) #unclear as the purpose of mod
print(dv)

if dv < pi:
	DM = 1
else:
	DM = -1
print(DM)

A = DM * sqrt(mag_r0 * mag_rf * (1 + cosDV))
print(A)

if dv == 0:
	print('Trajectory cant be computed.')

c2 = 0.5
c3 = .1666666666667
psi = 0
psi_up = 4 * exp(pi)
psi_low = -4 * pi
N = 0.8
count = 0
print(psi_low)

dt = dt0 + 1000
print(dt)

while (dt - dt0) > tol:
	y = mag_r0 + mag_rf + (A *(psi * c3 - 1))/sqrt(c2)

	if A > 0 and y < 0:
		psi = N/c3 *(1 - sqrt(c2)/A * (mag_r0 + mag_rf))

		if psi > 1e-6:
			c2 = (1 - cos(sqrt(psi)))/psi
			c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi*psi*psi)
		elif psi < -1e-6:
			c2 = (1 - cosh(sqrt(-psi)))/psi
			c3 = (sinh(sqrt(-psi)) - sin(sqrt(-psi)))/sqrt(-psi*-psi*-psi)
		else:
			c2 = 0.5
			c3 = .1666666666667

		y = mag_r0 + mag_rf + (A *(psi * c3 - 1))/sqrt(c2)
	print(mag_r0)
	print(mag_rf)
	print(psi)
	print(c2)
	print(c3)
	print(A)
	print(y)

	chi = sqrt(y/c2)
	dt = (chi^3 * c3 + A * sqrt(y))/sqrt(mu)

	if dt <= dt0:
		psi_low = psi
	else:
		psi_up = psi

	psi = 0.5 * (psi_up +psi_low)

	if psi > 1e-6:
		c2 = (1 - cos(sqrt(psi)))/psi
		c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi*psi*psi)
	elif psi < -1e-6:
		c2 = (1 - cos(sqrt(-psi)))/psi
		c3 = (sinh(sqrt(-psi)) - sin(sqrt(-psi)))/sqrt(-psi*-psi*-psi)
	else:
		c2 = 0.5
		c3 = .1666666666667

	eps = absolute(dt - dt0)

	count = count + 1

	if count > 1000:
		print('Lambert Solver hasn''t converged in 1000 iterations')
		break

f = 1 - y/mag_r0
gdot = 1 - y/mag_rf
g = A * sqrt(y/mu)

v0 = (rf - f * r0)/g
vf = (gdot * rf - r0)/g

print(v0)
print(vf)







































