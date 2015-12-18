import numpy
import math
import random

b = numpy.float_(0.2)

def psiT(r): #my trial wavefunction
    return (math.exp(-b*r))

def Elocal(r): 
    return -0.5*((b*b)-((2*b)/r))-(1/r)

def r(x,y,z): 
    return math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))

ne = 3
a = 0
e = numpy.array([random.random() for i in range (ne)])

dr = 1.05
"""
#equilibriate
for n in range (10000):
#move my electrons
    dt = numpy.array([dr*random.uniform(-1,1) for i in range (ne)])

    ep = e+dt
#was this move good?  
    R1 = r(e[0],e[1],e[2])
    Rp1 = r(ep[0],ep[1],ep[2])

    psip = psiT(Rp1)
    psi = psiT(R1)
    W = math.pow(psip/psi,2)  

    if W >= random.random():
        e = ep
        a = a+1
"""
#accumulate
El = 0
q = 0
srt = 0
for n in range (1000):
#move my electrons
    dt = numpy.array([dr*random.uniform(-1,1) for i in range (ne)])
    ep = e+dt
    #was this move good?  
    R1 = r(e[0],e[1],e[2])
    Rp1 = r(ep[0],ep[1],ep[2])
    print R1

    psip = psiT(Rp1)
    psi = psiT(R1)
    W = math.pow(psip/psi,2)    

    if W >= random.random():
        El = El+Elocal(Rp1)
        q = q+1
#ignore srt and s, I use them for error analysis which I have not     included here
        srt = srt+math.pow(Elocal(Rp1),2)
        s = math.pow((Elocal(Rp1))-(El/q),2)
        e = ep

print El/q
print q