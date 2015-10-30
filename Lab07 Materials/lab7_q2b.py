# PHY407, Fall 2015, Lab7, Q2b
# Author: DUONG, BANG CHI

# Import modules
from __future__ import division, print_function

from numpy import array, arange, empty, sqrt, float96, copy
from matplotlib.pyplot import figure, plot, scatter, annotate, axhline, \
     axvline, title, xlabel, ylabel, xlim, ylim, legend, show

# Define a class which solves ODEs using Bulirsch-Stoer method
class BulirschSolve:
	
    def __init__(self,f):
		
        self.f = f #self.array_decorator(f)

        self.initial_conditions = None
        self.solution = None
		
    def iterate(self, time_initial, time_final, H, delta):

        r0 = array(self.initial_conditions,float96)

        #H = (time_final - time_initial) / (365*10/7)

        tpoints = arange(time_initial, time_final, H)
        
        solution_rows = tpoints.shape[0]
        solution_cols = r0.shape[0]
        solution = empty([solution_rows, solution_cols], float)
		
	#r_points[0] = r0
        r = r0
        for i,t in enumerate(tpoints):

            solution[i] = r
            n = 1
            r1 = r + 0.5*H*f(r)
            r2 = r + H*f(r1)

            R1 = empty([n, solution_cols], float)
            R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

            error = 2*H*delta
            while error > H*delta:
                n += 1
                h = H/n

                r1 = r + 0.5*h*f(r)
                r2 = r + h*f(r1)
                for i in range(n-1):
                    r1 += h*f(r2)
                    r2 += h*f(r1)

                R2 = R1
                R1 = empty([n, solution_cols], float)
                R1[0] = 0.5 * (r1 + r2 + 0.5*h*f(r2))
                for m in range(1,n):
                    epsilon = (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
                    R1[m] = R1[m-1] + epsilon
                error = abs(epsilon[0])

            r = R1[n-1]

            self.H = H
            self.solution = solution
            self.t = tpoints

# Define the function f(r) = f(x,y,vx,vy)
def f(r):
	
	x,y,vx,vy = r
	
	Dx = vx
	Dy = vy
	
	R = sqrt(x**2 + y**2)
	
	Dvx = -G*M*x/R**3
	Dvy = -G*M*y/R**3
	
	return array([Dx,Dy,Dvx,Dvy],float96)

# Define constants
G = 6.6738e-11 # Newton's G in m^3/(kg.s^2)
M = 1.9891e30 # mass of the Sun in kg
m = 5.9722e24 # mass of the Earth in kg

# Define the problems that will be solved by Bulirsch-Stoer method
EarthOrbit = BulirschSolve(f)
PlutoOrbit = BulirschSolve(f)

# Initial conditions: [x0, y0, vx0, vy0]
EarthOrbit.initial_conditions = [1.4710e11, 0, 0, 3.0287e4]
PlutoOrbit.initial_conditions = [4.4368e12, 0, 0, 6.1218e3]

# Target accuracy delta: 1 km per year
delta_required = 1e3/365/24/60/60

# Solve it
EarthOrbit.iterate(0, 3600*24*365*10, H = 3600*24*7, delta = delta_required)
PlutoOrbit.iterate(0, 3600*24*365*1000, H = 3600*24*7e2, delta = delta_required)

# Call the solutions
x_EarthOrbit = EarthOrbit.solution[:,0]
y_EarthOrbit = EarthOrbit.solution[:,1]

x_PlutoOrbit = PlutoOrbit.solution[:,0]
y_PlutoOrbit = PlutoOrbit.solution[:,1]

# Plots
figure()
plot(x_EarthOrbit,y_EarthOrbit,'.', label='Eath\'s orbit')
plot(x_PlutoOrbit,y_PlutoOrbit,'o', label='Pluto\'s orbit')
legend().draggable()
sun = array([0, 0]) # Place sun at origin
scatter(sun[0], sun[1], 10, c=["r"])
xlim(-8e12, 8e12)
ylim(-8e12, 8e12)
axhline(color = 'black')
axvline(color = 'black')
xlabel('x (m)')
ylabel('y (m)')
title('Orbit of the Earth and Pluto around the Sun')

figure()
plot(x_EarthOrbit,y_EarthOrbit,'.', label='Eath\'s orbit')
xlim(-2e11, 2e11)
ylim(-2e11, 2e11)
scatter(sun[0], sun[1], 400, c=["r"])
annotate('Sun', xy = (0,0))
scatter(x_EarthOrbit[0], y_EarthOrbit[0], 50, c=["y"])
annotate('Perihelion', xy = (x_EarthOrbit[0]+0.25e10,y_EarthOrbit[0]))
axhline(color = 'black')
axvline(color = 'black')
xlabel('x (m)')
ylabel('y (m)')
title('Orbit of the Earth around the Sun in 10 years')

show()
