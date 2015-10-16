# PHY407, Fall 2015, Lab 4, Q2b
# Author: DUONG, BANG CHI

from scipy.constants import m_e, e, hbar
from numpy import tan, linspace

#constants
w = 1e-9 #m
V = 20*e #J (from eV)
accuracy = 0.001*e #eV

#define required functions
def tfun(E_in):
    return tan((w**2*m_e*E_in/(2*hbar**2))**0.5)

def EvenE(E_in):
    return ((V-E_in)/E_in)**0.5

def OddE(E_in):
    return - (E_in/(V-E_in))**0.5

def odd_fun(E_in):
    return tfun(E_in) - OddE(E_in)

def even_fun(E_in):
    return tfun(E_in) - EvenE(E_in)

#E in Joules for calculations
E = linspace(0,V, 1001)
#E in eV for plotting
E_ev = E/e

#------------------------------------------------------------------------------
# Define the binary search for Odd and Even Roots

def Odd_binary_search(a,b):
    while abs(b-a)>accuracy:
        midpoint = (a+b)/2
        #print midpoint
        if odd_fun(midpoint)*odd_fun(a) > 0:
            a = midpoint
        else:
            b = midpoint
    return midpoint/e

def Even_binary_search(a,b):
    while abs(b-a)>accuracy:
        midpoint = (a+b)/2
        #print midpoint
        if even_fun(midpoint)*even_fun(a) > 0:
            a = midpoint
        else:
            b = midpoint
    return midpoint/e

# The even roots
EvenRoot1 = Even_binary_search(0.3168*e, 0.3179*e)
EvenRoot2 = Even_binary_search(2.850*e, 2.853*e) 
EvenRoot3 = Even_binary_search(7.845*e, 7.855*e) 
EvenRoot4 = Even_binary_search(15.05*e, 15.10*e)
# The odd roots
OddRoot1 = Odd_binary_search(1.269*e, 1.271*e) 
OddRoot2 = Odd_binary_search(5.050*e, 5.052*e) 
OddRoot3 = Odd_binary_search(11.20*e, 11.23*e)
OddRoot4 = Odd_binary_search(19.12*e, 19.14*e)
print "Our first six roots are {}, {}, {}, {}, {}, {} eV".format(EvenRoot1,EvenRoot2,EvenRoot3,OddRoot1,OddRoot2,OddRoot3)

print "The Odd Roots (I found 4) are: {}, {}, {}, {}".format(OddRoot1, OddRoot2, OddRoot3, OddRoot4)
print "The Even Roots (I found 4) are: {}, {}, {}, {}".format(EvenRoot1, EvenRoot2, EvenRoot3, EvenRoot4)
