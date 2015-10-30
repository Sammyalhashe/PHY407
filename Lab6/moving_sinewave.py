from numpy import sin, linspace, pi
from pylab import plot, figure, axes, ylabel, xlabel, title, ion, pause

ion() # turn on interactive mode (so figure appears as code is running)

# grid
J=100
x=linspace(0,2*pi,J+1)

T=50 # number of timesteps
t=0

omega = 0.01 # angular frequency

y=sin(x-omega*t) # initial sine wave

# set up figure
fig = figure()
ax = axes(xlim=(0,2*pi), ylim=(-1,1))
#line, = ax.plot(x, y,'-b')
line = plot(x, y,'-b')
xlabel('x')
ylabel('y')
title('Timestep = '+str(t))


# loop through increasing t, recalculating and plotting sine wave at every timestep
for t in range(T):
    y=sin(x-omega*t)
    line[0].set_ydata(y)
    title('Timestep = '+str(t))
    pause(0.001)

