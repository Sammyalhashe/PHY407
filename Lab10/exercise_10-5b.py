# Use Mean Value Monte Carlo integration method to estimate integral of 1D
# function f(x) over specified range (a,b).
# Oliver Watt-Meyer, November 16, 2015. Adapted from Newman solutions.
from math import sin,sqrt,exp
from random import random

# define integrand
def f(x):
    return sin(1.0/(x*(2-x)))**2

N=10000 # number of points to use for Monte Carlo integration
sum_fx=0.0 # variable to keep track of sum(f(x))
sum_fxsq=0.0 # variable to keep track of sum(f(x)**2)

# domain of integration
a=0.0
b=2.0

for i in range(N):
    # pick random value between a and b
    x=(b-a)*random()+a
    fx=f(x)
    
    sum_fx+=fx
    sum_fxsq+=fx**2
    
# compute mean and variance of series of values fx
mean=sum_fx/N
var=sum_fxsq/N - mean**2

# Compute integral and error using Eqs. 10.30 and 10.32
I=(b-a)*mean
sigma=(b-a)*sqrt(var/N)

print 'With Mean Value method, I = %.4f plus/minus %.4f'%(I,sigma)