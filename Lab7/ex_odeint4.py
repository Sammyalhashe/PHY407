from pylab import *
from numpy import *

#RHS vector of ODE
def f(y,t):
    return -2.3*y

#Analytic solution
def analytic(t,c):
    return c*exp(-2.3*t)


#t for analytic solution
endt=5
t=linspace(0.0,float(endt),20)

#initial condition
y=zeros(20)
y0=1.0
#calculate analytic solution
y=analytic(t,y0)


#integrate ODE using forward Euler method for unstable h
h_un=1.0
N_un=endt//h_un+1
y_un=zeros(N_un)
y_un[0]=1.0
t_un =zeros(N_un)
t_un[0]=0.0

for i in range(len(y_un)-1):
    t_un[i+1]=h_un+t_un[i]
    y_un[i+1]=y_un[i]+h_un*f(y_un[i],t_un[i])

#integrate ODE using forward Euler method for stable h
h_st=0.7
N_st=endt//h_st+1
y_st=zeros(N_st)
y_st[0]=1.0
t_st =zeros(N_st)
t_st[0]=0.0

for i in range(len(y_st)-1):
    t_st[i+1]=h_st+t_st[i]
    y_st[i+1]=y_st[i]+h_st*f(y_st[i],t_st[i])


#plot analytic solution
plot(t,y,'-',label="analytic")


#plot unstable odeint solution
plot(t_un,y_un,'o-',label="unstable h=1")

#plot stable odeint solution
plot(t_st,y_st,'s-',label="stable h=0.7")

title("dy/dt=-2.3y with forward Euler")
legend(loc=2)
savefig("fig_ex_odeint4.jpg")

show()
