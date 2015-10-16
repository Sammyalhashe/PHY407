# PHY407, Fall 2015, Lab 4, Q1c
# Author: DUONG, BANG CHI

# Import modules
from numpy import linspace, pi, exp, array
from matplotlib.pyplot import figure, plot, legend, title, xlabel, ylabel, show
from numpy.linalg import solve

# Define constants
# Resistance in ohm
R1, R3, R5 = 1e3, 1e3, 1e3
R2, R4, R6 = 2e3, 2e3, 2e3

# Capacitance in F
C1 = 1e-6
C2 = 0.5e-6

x_positive = 3.0 # in V

w = 1e3 # in 1/s

# The matrix for the three equation
A = array([[ 1/R1 + 1/R4 + complex(0, w*C1), complex(0, -w*C1)                  , 0                              ],
           [ complex(0, -w*C1)             , 1/R2 + 1/R5 + complex(0, w*(C1+C2)), complex(0, -w*C2)              ],
           [ 0                             , complex(0, -w*C2)                  , 1/R3 + 1/R6 + complex(0, w*C2) ]])
v = array([x_positive/R1, x_positive/R2, x_positive/R3])
x_solution = solve(A,v)
print x_solution

# Print out V1, V2, V3 based on x_solution
from cmath import polar
for i,xi in enumerate(x_solution):
    r,theta = polar(xi)
    print("{}: V{} = {}, phase = {}".format(i, i+1, r, theta*180/pi))

# Define time
time = linspace(0, 4.0*pi/w, 101) # a cycle through two periods

# Define a function V(x,t)
def V(x_value, time_value):
    return x_value * exp(1j*w*time_value)

# Plots
for k in range(len(x_solution)):
    plot(time, V(x_solution[k], time).real, label = "V%s" %(k+1))
legend().draggable()
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage as a function of time')

show()
