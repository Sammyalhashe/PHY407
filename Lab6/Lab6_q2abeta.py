#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

"""
Using RK4, solves and plots the one dimensional 
system of the N identical masses with springs
"""

__author__ = "Eric Yeung"


m = 1.
k = 6.
omega = 2.
N = 5

# Define the driving force

def F_i(i, t):

	if i == 1:
		drivingforce = np.cos(omega*t)

	else:
		drivingforce = 0

	return drivingforce


def f(r, t):
	
	#fArray = np.empty([2, N], float)

	for i in range(N):

                zeta = r[0]
                v = r[1]

		# Remember to shift indices
                
		dzeta = v

		if i == 1:
			dvi = k/m*(zeta[2] - zeta[1]) + F_i(1, t)

		elif i == N:
			dvN = k/m*(zeta[N-1] - zeta[N]) # F_N = 0

		else:
			dv[i] = k/m*(zeta[i+1] - zeta[i] + zeta[i-1] - zeta[i]) # F_i = 0

		#temp = fList.append(dzeta, dv[1], dv[i], dv[N]) # Should be length 2*N
		#fArray = np.array(temp)

		#fArray[i, :] = dsi # ith row

	return np.array([dzeta, dv[1], dv[i], dv[N]], float)


tplot = np.arange(0, 21)



"""
def RK4(h, s, f):
	k1 = f(t,u)
	k2 = f(t+dt*0.5,u+k1*0.5*dt)
	k3 = f(t+dt*0.5,u+k2*0.5*dt)
	k4 = f(t+dt,u+k3*dt)



	v += dt * (k1+2*k2+2*k3+k4)/6

	# v doesn't explicitly depend on other variables
	k1 = k2 = k3 = k4 = v

	u += dt * (k1+2*k2+2*k3+k4)/6

	return u,v

"""

"""

		def zeta(i):

			zeta1 = r[0] # i = 1
			zeta2 = r[1] # i = 2, ...

			return

		def s(i):

			return



for j in range(N):
	figure(j)
	plt.plot



u' = v
v' = f(t,u)

def rk_iter(u, v, t, dt):
    k1 = f(t,u)
    k2 = f(t+dt*0.5,u+k1*0.5*dt)
    k3 = f(t+dt*0.5,u+k2*0.5*dt)
    k4 = f(t+dt,u+k3*dt)

    v += dt * (k1+2*k2+2*k3+k4)/6

    # v doesn't explicitly depend on other variables
    k1 = k2 = k3 = k4 = v

    u += dt * (k1+2*k2+2*k3+k4)/6

    return u,v




"""
