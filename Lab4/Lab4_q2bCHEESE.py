#!/usr/bin/env python
from __future__ import division
from Lab4_q2b_FUNCTIONS import *

"""
Finds the first six energy levels using a binary search
"""

__author__ = "Eric Yeung"


tol = 1000000
range1 = np.linspace(0.317*e, 0.318*e, tol)

E1 = []
F1 = []

for i in range(tol):
	E1.append(EvenE(V, range1[i]))
	F1.append(tfun(range1[i], m_e, w, hbar))

def intersect(a, b):
    result=[]
    for i in b:
        if isinstance(i,list):
            result.append(intersect(a,i))
        else:
            if i in a:
                 result.append(i)
    return result

print intersect(E1, F1)