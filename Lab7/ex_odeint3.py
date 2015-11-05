from scipy.integrate import odeint
from pylab import *

#RHS vector of ODE
def f(y,t):
    return -1.2*y

#Analytic solution
def analytic(t,c):
    return c*exp(-1.2*t)


#t for analytic solution
t=linspace(0.0, 2.0,20)

#initial condition
y0=linspace(0.5, 1.5, 5)
y=zeros([len(y0),20])

#calculate analytic solution
for i in range(len(y0)):
    y[i,:]=analytic(t,y0[i])

#time array 
ts=[0.0,0.5,1.0,1.5, 2.0]

#integrate ODE
ys=odeint(f,1.,ts)

#plot analytic solution
for i in range(len(y)):
    plot(t,y[i],'-')


#plot odeint solution
plot(ts,ys,'o-',label="odeint")

title("dy/dt=-1.2y")
savefig("fig_ex_odeint3.jpg")
show()
