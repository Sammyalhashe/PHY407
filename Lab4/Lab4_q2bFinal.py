#!/usr/bin/env python
from __future__ import division
from Lab4_q2b_FUNCTIONS import *

"""
Finds the first six energy levels using a binary search
"""

__author__ = "Eric Yeung"

tol = 1.60217662e-22 # accuracy of 0.001 eV

def binarysearch_EVEN(ELower, EUpper):
	FLower = EvenE(V, ELower)
	EMid = (ELower + EUpper)/2.
	
	FMid = EvenE(V, EMid)
	Intersect = tfun(EMid,m_e,w,hbar)

	i = 0
	while abs(EUpper - ELower) >= tol:
		
		i += 1 # Count the iterations
		if FMid > Intersect:
			EUpper = EMid

		elif FMid < Intersect:
			ELower = EMid

		EMid = (ELower + EUpper)/2.
		FMid = EvenE(V, EMid)
		Intersect = tfun(EMid,m_e,w,hbar)

	return Intersect, FMid, EMid, i, EUpper, ELower

def binarysearch_ODD(ELower, EUpper):
	FLower = OddE(V, ELower)
	EMid = (ELower + EUpper)/2.
	
	FMid = OddE(V, EMid)
	Intersect = tfun(EMid,m_e,w,hbar)

	i = 0
	while abs(EUpper - ELower) >= tol:
		
		i += 1 # Count the iterations
		if FMid > Intersect:
			EUpper = EMid

		elif FMid < Intersect:
			ELower = EMid

		EMid = (ELower + EUpper)/2.
		FMid = OddE(V, EMid)
		Intersect = tfun(EMid,m_e,w,hbar)

	return Intersect, FMid, EMid, i, EUpper, ELower

if __name__ == "__main__":

    #Starting with the even roots

    FirstRoot = binarysearch_EVEN(0.317*e, 0.318*e) # First root
    SecondRoot = binarysearch_EVEN(0.3625*e, 0.3630*e) # Second root, this value is too close to asymptote
    FifthRoot = binarysearch_EVEN(2.850*e, 2.855*e) # Fifth root
    #print FirstRoot; print SecondRoot; print FifthRoot

    # The odd roots

    ThirdRoot = binarysearch_ODD(0.3638*e, 0.3640*e) # Third root
    FourthRoot = binarysearch_ODD((0.00004 + 1.27008)*e, (0.00005 + 1.27008)*e) # Fourth root, dont know why this is wrong.
    SixthRoot = binarysearch_ODD(5.050*e, 5.055*e) # Sixth root
    #print ThirdRoot; FourthRoot; print SixthRoot

    print "Our first six energy levels are %s, %s, %s, %s, %s, %s eV" % (FirstRoot[2]/e, SecondRoot[2]/e, ThirdRoot[2]/e, FourthRoot[2]/e, FifthRoot[2]/e, SixthRoot[2]/e)
