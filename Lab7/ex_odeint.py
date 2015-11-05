from scipy.integrate import odeint
from pylab import *

#RHS vector of ODE
def f(y,t):
    return 1.2*y

#Analytic solution
def analytic(t,c):
    return c*exp(1.2*t)


#t for analytic solution
t=linspace(0.0,2.0,20)

#initial condition
y0=1.0

#calculate analytic solution
y=analytic(t,y0)

#time array 
ts=[0.0,0.5,1.0, 1.5,2.0]
ys=odeint(f,1.,ts)

#plot analytic solution
plot(t,y,'-',label="analytic")
#plot odeint solution
plot(ts,ys,'o-',label="Euler")

legend(loc="upper left")
title("dy/dt=1.2y")
savefig("fig_ex_odeint.jpg")
show()
