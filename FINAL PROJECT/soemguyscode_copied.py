import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1337) # Seed the random generator 

def psiT(b, r): #my trial wavefunction
    return (np.exp(-b*r))

def Elocal(b, r): 
    return -0.5*((b*b)-((2*b)/r))-(1/r)

def r(x,y,z): 
    return np.sqrt(x**2+y**2+z**2)

def main(b):
    ne = 3
    a = 0
    e = np.array([np.random.rand() for i in range (ne)])

    dr = 1.05
    #accumulate
    El = 0
    q = 0
    srt = 0
    for n in range (1000):
    #move my electrons
        dt = np.array([dr*np.random.uniform(-1,1) for i in range (ne)])

        ep = e+dt
        #was this move good?  
        R1 = r(e[0],e[1],e[2])
        Rp1 = r(ep[0],ep[1],ep[2])

        psip = psiT(b, Rp1)
        psi = psiT(b, R1)
        W = (psip/psi)**2    

        if W >= np.random.rand():
            El = El+Elocal(b, Rp1)
            q = q+1
            e = ep

    return El/q

brange = np.arange(0.2, 5.5, 5.3/1000)
# Initialize local energy and variance lists
energylist = []
# Calculate the ground state for varying alpha
for b in brange:
    energylist.append(main(b))

plt.plot(brange, energylist, color = 'r', linestyle = '--', label = 'VMC Simulation')
plt.axhline(y = -0.5, label = "Exact solution")
plt.ylabel('Ground State Energy')
plt.legend().draggable()
plt.show()