# PHY407, Fall 2015, Lab09, Q2d
# Author: DUONG, BANG CHI

# Import modules
from numpy import array, empty, diag, exp, arange, dot, real,copy,zeros
from numpy.linalg import solve
from scipy.constants import hbar
from matplotlib.pyplot import *

# Define constants
m = 9.109e-31 # electron mass in kg
L = 1e-8 # box length in m
N = 1000 # spatial slices
a = L/N # distance between grid points
h = 1e-18 # size of time step in seconds
t_end = 1e-15

# Define wave function at t=0
def psi0(x):
    x0 = 0.4*L # in m
    sigma = 1e-10 # in m
    k = 5*10**(10) # in m^(-1)
    return exp(-((x-x0)**2.0)/(2.0*sigma**2))*exp(1j*k*x)

# Define V(x)
def V(x):
    if x <= L/2:
        return 0.0
    elif x > L/2:
        return 40.0 * 1.60218e-19

x_array = arange(0,L,a) # x points from 0 to L


# Tridiagonal Matrices A and B
a1 = zeros([N],complex)
b1 = zeros([N],complex)

for i in range(len(x_array)):
    a1[i] = (1 + 1j*h*(hbar)/(2.0*m*a**2) + V(x_array[i])*1j/(2*hbar)) # components of A
    b1[i] = (1 - h*(1j*hbar)/(2.0*m*a**2) - V(x_array[i])*1j/(2*hbar)) # components of B
a2 = -h*(1j*hbar/(4.0*m*a**2))
b2 = h*(1j*hbar/(4.0*m*a**2))

def tridiag(a0,b0,c0):
    return diag(a0,-1)+diag(b0,0)+diag(c0,1)

A = tridiag([a2]*(N-1),a1,[a2]*(N-1)) # Matrix A
B = tridiag([b2]*(N-1),b1,[b2]*(N-1)) # Matrix B

t_array = arange(0,t_end,h) # time array in 1 seconds with timestep h

v0 = dot(B, psi0(x_array))
psi1 = solve(A,v0)


# Define psi
psi = []
for t in t_array:
    v1 = dot(B,psi1)
    psi2 = solve(A, v1)
    psi1= copy(psi2)
    psi.append(psi1)
psi = array(psi,complex)

# The animation
ion() #starts animation
fig = figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x_array,abs(psi[0]),label = '$|\Psi(x,t)|$')
line2, = ax.plot(x_array,psi[0].real,label = '$\Psi(x,t)$')
line3, = ax.plot(x_array,-abs(psi[0]),label = '$-|\Psi(x,t)|$')
xlabel('$x$')
ylabel('$\Psi(x,t)$')
legend()
draw()
t = 0
count = 0
while t < t_end:
    title('when t = '+str(t))
    line1.set_ydata(abs(psi[count]))
    line2.set_ydata(psi[count].real)
    line3.set_ydata(-abs(psi[count]))
    draw()
    pause(0.01)
    t += h
    count += 1
ioff() # ends animation
show()
