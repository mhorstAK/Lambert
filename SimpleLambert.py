#!/usr/bin/env python
import math
import numpy as np
from math import atan2, pi, sqrt, sin, cos, cosh
from numpy import absolute, dot, array, sqrt, linalg, mod


r0 = [-96948447.3751163, 46106976.1901101, 6225649.83906096]
rf = [31161153.5750633, 143995536.213184, 3018.4216823707]
dt0 = 26179200

mu = 132712440018 #gravitational parameter of the sun
tol = .000001 #tolerance thingy (I think)

psi = 0
c2 = 0.5 
c3 = .16666666666667 
psi_up  = 4*(pi**2) 
psi_low = -4*pi
N = 0.8

mag_r0 = []
mag_rf = []
dnu = []
DM = 0
cdnu = 0



def cdnuYup(r0, rf):
	mag_r0 = np.linalg.norm(r0)
	mag_rf = np.linalg.norm(rf)
	cdnu = dot(r0,rf)/(mag_r0*mag_rf)
	return mag_r0
	return mag_rf
	return cdnu

def moDD(r0, rf):
	nu1 = atan2(r0[1],r0[0])
	nu2 = atan2(rf[1],rf[0])
	dnu = nu2 - nu1
	dnu = mod(dnu,2*pi)
	return dnu

def direction(dnu, pi):
	if dnu < pi:
		DM = 1
	else:
		DM = -1
	return DM

def AAA(DM, mag_r0, mag_rf, cdnu):
	A = DM*sqrt(mag_r0*mag_rf*(1+cdnu))
	return A

def dnuCheck(dnu):
	if dnu==0: 
		print('Delta Nu = 0!!!')

def psiTol(psi, tol):
	if psi > tol:
		c2 = (1 - cos(sqrt(psi)))/psi
		c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi**3)
	elif psi < -tol:
		c2 = (1 - cosh(sqrt(-psi)))/psi
		c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)**3)
	else:
		c2 = 0.5
		c3 = .16666666666667
	return c2
	return c3

def AandY(N, c2, c3, A, y, tol, mag_r0, mag_rf):
	if A > 0 and y < 0:
		psi = N/c3 * (1 - sqrt(c2)/A * (mag_r0+mag_rf))
		psiTol(psi, tol)
		y = mag_r0 + mag_rf + (A*(psi*c3 - 1))/sqrt(c2)
	return y
	return psi
	return c2
	return c3

def dtDt0(dt, dt0, psi):
	if dt <= dt0:
		psi_low = psi
	else:
		psi_up = psi
	return psi_up
	return psi_low

def theWhile(eps, N, c2, c3, A, tol, mag_r0, mag_rf, mu, dt0, while_count):
	while eps > tol:
		y = mag_r0 + mag_rf + (A*(psi*c3 - 1))/sqrt(c2)
		AandY(N, c2, c3, A, y, tol, mag_r0, mag_rf)
		chi = sqrt(y/c2)
		dt = ((chi**3)*c3 + A*sqrt(y))/sqrt(mu)
		dtDt0(dt, dt0, psi)
		psi = 0.5*(psi_up + psi_low)
		psiTol(psi, tol)
		eps = absolute(dt - dt0)
		while_count = while_count + 1
		if while_count > 1000:
			break



def main():
	cdnuYup(r0, rf)
	moDD(r0, rf)
	direction(dnu, pi)
	AAA(DM, mag_r0, mag_rf, cdnu)
	dnuCheck(dnu)
	eps = dt0 + 1000
	while_count = 0
	theWhile(eps, N, c2, c3, A, tol, mag_r0, mag_rf, mu, dt0, while_count)

	f = 1 - y/mag_r0
	g_dot = 1 - y/mag_rf
	g = A*sqrt(y/mu)

	v0 = (np.array(rf) - np.multiply(f,r0))/g
	vf = (np.multiply(g_dot,rf) - r0)/g
	print(v0)
	print(vf)

main()