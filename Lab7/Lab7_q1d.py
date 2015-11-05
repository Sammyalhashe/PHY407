#!/usr/bin/env python
from __future__ import division
from Lab7_q1ab import *
from Lab7_q1c import *

"""
Compares computational and Analytical
"""

__author__ = "Eric Yeung"

# For l = 0, n = 1:
def R_01(r):
    return 2/(a_0**(3/2))*np.exp(-r/a_0)

def R_02(r):
    return 1/(2*np.sqrt(2)*a_0**(3/2))*(2 - r/a_0)*np.exp(-r/(2*a_0))

# For l = 1, n = 1, we have degeneracy. For l = 1, n = 2:
def R_12(r):
    return 1/(2*np.sqrt(6)*a_0**(3/2))*(r/a_0)*np.exp(-r/(2*a_0))

# ANALYTICAL
plt.plot(rplot, R_01(rplot), label = 'Ground state l = 0')
plt.plot(rplot, R_02(rplot), label = '1st Excited State l = 0')
plt.plot(rplot, R_12(rplot), label = 'Ground state l = 1')

# COMPUTED
plt.plot(rplot, norm01, label = 'Ground state l = 0')
plt.plot(rplot, norm02, label = '1st Excited State l = 0')
plt.plot(rplot, norm11, label = 'Ground state l = 1')

plt.xlabel('r')
plt.ylabel('R(r)')
plt.title('Computed Wavefunctions vs Analytical Solutions')
plt.legend().draggable()
	
plt.show()
