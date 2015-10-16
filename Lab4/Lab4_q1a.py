#!/usr/bin/env python
from __future__ import division

#Python module for PHY407
#Modified the script by Paul Kushner, 2015-09-26

import numpy as np
#The following will be useful for partial pivoting
#from numpy import empty, copy

#Implement Gaussian Elimination
#This should be non-destructive for input arrays, so we will copy A and v to temporary variables

def GaussElim(A_in,v_in):
    #copy A and v to temporary variables using copy command
    A = np.copy(A_in)
    v = np.copy(v_in)
    N = len(v)
    
    for m in range(N):
    
        # Divide by the diagonal element
        div = A[m,m]
        A[m,:] /= div
        v[m] /= div
    
        # Now subtract from the lower rows
        for i in range(m+1,N):
            mult = A[i,m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]
        
    # Backsubstitution
    #create an array of the same type as the input array
    x = np.empty(N,dtype=v.dtype)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m,i]*x[i] 
    return x

def PartialPivot(A_in,v_in):
    A = np.copy(A_in)
    v = np.copy(v_in)
    N = len(v)
    for i in range(N):

        div = A[i,i]
        A[i,:] /= div
        v[i] /= div

        for j in range(i+1,N):
            
            # Swap rows so we get the greatest value
            if abs(A[i,j]) > abs(A[j,i]):
                A[j,:], A[i,:] = np.copy(A[i,:]), np.copy(A[j,:])
            #A[j,:], A[i,:] = np.copy(A[i,:]), np.copy(A[j,:])   
                v[j], v[i] = v[i], v[j]
                print A

            #mult = A[j,i]
            #A[j,:] -= mult*A[i,:]
            #v[j] -= mult*v[i]
            
    # Backsubstitution
    #create an array of the same type as the input array
    x = np.empty(N,dtype=v.dtype)
    for i in range(N-1,-1,-1):
        x[i] = v[i]
        for i in range(i+1,N):
            x[i] -= A[i,j]*x[j] 
    return x

if __name__ == '__main__':
    
    A = np.array([[ 2,  1,  4,  1 ],
                  [ 3,  4, -1, -1 ],
                  [ 1, -4,  1,  5 ],
                  [ 2, -2,  1,  3 ]], float)
    v = np.array([ -4, 3, 9, 7 ],float)

    print "The known answer from 6.16 is"
    print GaussElim(A, v)

    print "The partial pivoting method gives"
    print PartialPivot(A,v)    

    print "They agree!"
