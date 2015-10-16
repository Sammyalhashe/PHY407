#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from gaussxw import gaussxw

"""
Temperature of light bulb
"""

__author__ = "Eric Yeung"

def efficiency(T):

	N = 100
	planck = 6.62607004e-34 # Joules * s
	boltz = 1.38064852e-23 # Joules / K
	c = 299792458 # m / s
	
	lambda1 = 3.9e-7 # metres
	lambda2 = 7.5e-7 # metres
	x0 = planck*c/(lambda2*boltz*T)
	xf = planck*c/(lambda1*boltz*T)

	def I(x):
		return (x**3)/(np.exp(x) - 1) # overflow?

	x,w = gaussxw(N)
	xp = 0.5*(xf-x0)*x + 0.5*(xf+x0)
	wp = 0.5*(xf-x0)*w

	s = 0.0
	for k in range(N):
		s += wp[k]*I(xp[k])
	
	return 15/(np.pi**4)*s


if __name__ == '__main__':
	
	Tplot = np.linspace(300, 10000, 100)
	EffArray = []

	for i in range(100):
		EffArray.append(efficiency(Tplot[i]))

	plt.plot(Tplot, EffArray)
	plt.xlabel('Temperature (K)')
	plt.ylabel('Efficiency')
	plt.title('Efficiency for T = [300, 10000]')

	plt.show()


"""
Finding the maximum using the golden ratio search
"""

Tolerance = 1e-4 # Accuracy of 1K
gold = (1+np.sqrt(5))/2

t1 = 6000.
t4 = 8000.
t2 = t4 - (t4 - t1)/gold
t3 = t1 + (t4 - t1)/gold

n1 = -1*efficiency(t1)
n2 = -1*efficiency(t2)
n3 = -1*efficiency(t3)
n4 = -1*efficiency(t4)

while t4 - t1 > Tolerance:
	if n2 < n3:
		t4, n4 = t3, n3
		t3, n3 = t2, n2
		t2 = t4 - (t4 - t1)/gold
		n2 = -1*efficiency(t2)

	else:    
		t1, n1 = t2, n2
		t2, n2 = t3, n3
		t3 = t1 + (t4 - t1)/gold
		n3 = -1*efficiency(t3)

print "The maximum efficiency %.5f occurs at %.2f K" % (efficiency(0.5*(t1 + t4)), 0.5*(t1 + t4))
