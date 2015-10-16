# PHY407, Fall 2015, Lab 4, Q2a
# Author: DUONG, BANG CHI

from numpy import tanh, cosh, linspace
from pylab import figure, subplot, plot, show, title, ylim, xlabel, ylabel, legend
import scipy.optimize

Tmax = 2.0
points = 1000
accuracy = 1e-6

mag_relaxation = []
mag_newton = []
iter_relaxation = []
iter_newton = []
temp = linspace (0.01,Tmax,points)

#-------------------------Relaxation method-----------------------
for T in temp:
    m1 = 1.0
    error = 1.0
    iter_num = 0
    while error>accuracy:
        m1,m2 = tanh(m1/T),m1
        error = abs((m1-m2)/(1-T*cosh(m1/T)**2))
        iter_num += 1

    mag_relaxation.append(m1)
    iter_relaxation.append(iter_num)


#---------------------------Newton's method------------------------

for T in temp:
    m = 1.0
    delta = 1.0
    iter_num = 0
    while abs(delta)>accuracy:
        delta = (m - tanh(m/T))/(1/T*cosh(m/T)**2)
        m -= delta
        iter_num += 1
    
    iter_newton.append(iter_num)
    mag_newton.append(m)

	
# Plot Magnetization vs Temperature
figure()
plot(temp, mag_relaxation, label='Relaxation method')
plot(temp, mag_newton, label='Newton method')
ylim(-0.1, 1.1)
xlabel('Temperature')
ylabel('Magnetization')
legend(loc='upper right')


# Plot Number of Iteration for 2 methods
figure()
subplot(211)
plot(temp, iter_newton)
title("Number of iterations for Newton's method")
ylabel("Count")

subplot(212)
plot(temp, iter_relaxation)
title("Number of iterations for Relaxation method")
xlabel("Temperature")
ylabel("Count")

#------------------------------------------------------------------
# Verify the power law with beta ~ 0.5
def powerlaw(T, slope, const, beta):
    return slope*(1-T)**beta + const

# Make a new array for temperature that is very near 1.0
new_temp = linspace(0.97,1.0,points)
new_mag_newton=[]
# Re-define newton's method for the new temperature array
for T in new_temp:
    m = 1.0
    delta = 1.0
    iter_num = 0
    while abs(delta)>accuracy:
        delta = (m - tanh(m/T))/(1/T*cosh(m/T)**2)
        m -= delta
    new_mag_newton.append(m)
# Fitted line with beta(guess) = 0.5
popt,pcov = scipy.optimize.curve_fit(powerlaw,new_temp,new_mag_newton,p0=(1,0,0.5))

powerlaw_curve = powerlaw(new_temp,popt[0],popt[1],popt[2])
# Plot and compare
figure()
plot(new_temp, new_mag_newton, label='Newton method')
plot(new_temp, powerlaw_curve,label='Powerlaw function with beta = {}'.format(popt[2]))
ylim(-0.1, 0.4)
xlabel('Temperature')
ylabel('Magnetization')
legend(loc='upper right')

show()


