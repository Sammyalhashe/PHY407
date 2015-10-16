#!/usr/bin/env python
from __future__ import division
from Lab4_q2b_a import *

"""
Finds the first sex energy levels using a binary seartch
"""

__author__ = "Eric Yeung"


Tolerance = 1e-54
E0 = 0.31*e # Joules
E1 = 0.32*e # Joules

mid = 0.5*(E0 + E1)
tanmid = tfun(mid,m_e,w,hbar)
fmid = EvenE(V, mid)
#print fmid, tanmid


#print 5.04685617975e-20/e, EvenE(V, 5.04685617975e-20), tfun(4.9667e-20,m_e,w,hbar)

while abs(E1 - E0) > Tolerance:

	if EvenE(V, mid) == tfun(mid,m_e,w,hbar):
		print mid

	elif EvenE(V, mid) < tfun(mid,m_e,w,hbar):
		E0 = mid
		print "lower"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)

	elif EvenE(V, mid) > tfun(mid,m_e,w,hbar):
		E1 = mid
		print "higher"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)
		print EvenE(V, mid), tfun(mid,m_e,w,hbar)

	else:
		pass


"""

def binarysearch_EVEN(ELower, EUpper, tol):
	EMid = (ELower + EUpper)/2.
	
	FMid = EvenE(V/e, EMid)
	Intersect = tfun(EMid,m_e,w,hbar)

	i = 0

	while abs(EUpper - ELower) > tol:
		
		i += 1

		if FMid > Intersect:
			ELower = EMid

		elif FMid < Intersect:
			EUpper = EMid

		else:
			pass

		EMid = (ELower + EUpper)/2.
		FMid = EvenE(V/e, EMid)

	return tfun(EMid,m_e,w,hbar), FMid, EMid, i, EUpper, ELower

print binarysearch_EVEN(0.31, 0.32, 1e-444)


"""


"""

while EvenE(V, mid) != tfun(mid,m_e,w,hbar):

	if fmid < tfun(mid,m_e,w,hbar):
		E0 = mid
		print "less than"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)

	elif fmid > tfun(mid,m_e,w,hbar):
		E1 = mid
		print "more than"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)

	#print EvenE(V, mid)


while abs(E1 - E0) > Tolerance:

	mid = 0.5*(E0 + E1)

	if abs(EvenE(V, mid)/tfun(mid,m_e,w,hbar)) < Tolerance: 
		print mid

	elif EvenE(V, mid) < tfun(mid,m_e,w,hbar):
		E0 = mid
		print "less than"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)

	elif EvenE(V, mid) > tfun(mid,m_e,w,hbar):
		E1 = mid
		print "mid"
		mid = 0.5*(E0 + E1)
		EvenE(V, mid)

	else:
		pass

	print EvenE(V, mid)

        
def binary_search(seq, t):
    min = 0
    max = len(seq) - 1
    while True:
        if max < min:
            return -1
        m = (min + max) // 2
        if seq[m] < t:
            min = m + 1
        elif seq[m] > t:
            max = m - 1
        else:
            return m

if EvenE(V, mid) == tfun(mid,m_e,w,hbar):

"""
