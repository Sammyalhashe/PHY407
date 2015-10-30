from pylab import *
# set constants

#dv/dt = -(k/m) x
#dx/dt = v
#v(t + dt) = v(t) - dt (k/m) x(t)
#x(t + dt) = x(t) + dt v(t)

k = 10
m = 1
t_max = 10.0
#Set the number of iterations to be used in the for loop
no_of_iterations=1000
 
# set time step so that the loop will always iterate until t=t_max seconds 
dt = t_max/no_of_iterations
 
 
# make arrays to store data
t = zeros(no_of_iterations)
x = zeros(no_of_iterations)
v = zeros(no_of_iterations)
 
# set initial conditions
t[0] = 0
x[0] = 5
v[0] = 0
 
# evolve
for i in arange(1,no_of_iterations):
    t[i] = dt * i
    v[i] = v[i-1] - dt *k/m*x[i-1]
    x[i] = x[i-1] + dt * v[i-1]
    

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(t, x)
ax1.set_title('Position Graph')
ax2.plot(t, v, color='r')
ax2.set_title('Velocity Graph')
ax1.set_xlabel('Time')
ax2.set_xlabel('Time')
ax1.set_ylabel('Position')
ax2.set_ylabel('Velocity')

plt.show()