#!/usr/bin/env python
from Lab6_q3c import *
"""
Places dots showing the positions
Should be more closely spaced near the perihelion
And further spaced near the aphelion
"""

__author__ = "Eric Yeung"


# Plotting x vs y
plt.plot(r_sols[0,:], r_sols[1,:], '.')
plt.annotate('SUN', xy=(-1,0), bbox=dict(boxstyle="circle", fc="r"))
plt.xlabel('x')
plt.ylabel('y')
plt.title('Adaptive Step Size: Positions')
plt.show()