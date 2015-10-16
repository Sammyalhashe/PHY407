#PHY407 2015
#Solution to Lab 04, Q2a
#Paul Kushner
from scipy.constants import m_e, e, hbar
import numpy as np
import matplotlib.pyplot as plt

#define required functions
def tfun(E_in,m_e,w,hbar):
    return np.tan((w**2*m_e*E_in/(2*hbar**2))**0.5)

def EvenE(V,E_in):
    return ((V-E_in)/E_in)**0.5

def OddE(V,E_in):
    return - (E_in/(V-E_in))**0.5

#define constants
w = 1e-9 #m
V = 20*e #J (from eV)

#E in Joules for calculations
E = np.linspace(V*1e-6,V, 1001)
#E in eV for plotting
E_ev = E/e

#plot the functions over four different ranges, to help identify initial E's for the bisection solver.

if __name__ == "__main__":
    
    plt.figure(1)
    plt.plot(E/e, tfun(E, m_e, w, hbar), E/e, EvenE(V,E), E/e, OddE(V,E))
    plt.legend(('tan function', 'Even', 'Odd'))
    plt.xlabel('E(eV)')
    plt.ylim((-10,10))

    plt.figure(2)
    plt.plot(E/e, tfun(E, m_e, w, hbar), E/e, EvenE(V,E), E/e, OddE(V,E))
    plt.legend(('tan function', 'Even', 'Odd'),)
    plt.xlabel('E(eV)')
    plt.xlim((0.2,0.5))
    plt.ylim((5,10))

    plt.figure(3)
    plt.plot(E/e, tfun(E, m_e, w, hbar), E/e, EvenE(V,E), E/e, OddE(V,E))
    plt.legend(('tan function', 'Even', 'Odd'))
    plt.xlabel('E(eV)')
    plt.xlim((0.2,0.5))
    plt.ylim((-2,0.0))

    plt.figure(4)
    plt.plot(E/e, tfun(E, m_e, w, hbar), E/e, EvenE(V,E), E/e, OddE(V,E))
    plt.legend(('tan function', 'Even', 'Odd'))
    plt.xlabel('E(eV)')
    plt.xlim((1,3))
    plt.ylim((-2,4))

    plt.figure(5)
    plt.plot(E/e, tfun(E, m_e, w, hbar), E/e, EvenE(V,E), E/e, OddE(V,E))
    plt.legend(('tan function', 'Even', 'Odd'))
    plt.xlabel('E(eV)')
    plt.xlim((4.5,6))
    plt.ylim((-2,0))

    plt.show()

else:
    print "The plots are being omitted!"

