#quick demo of leapfrog method for reference
#using linear spring F=-kx=ma
#we could solve this using Verlet and save time but I used leapfrog to give you a working example.
from pylab import plot, show
#initialize $\eta$ 
def f(x,v):
   return v,-k*x/m
nsteps = 200
h=0.3
k=1.0
m=1.0
t=0
x_full = 1.0
v_full = 0.0
plot(t, x_full,'b+')
dx,dv = f(x_full, v_full)
x_half = x_full + 0.5*h*dx
v_half = v_full + 0.5*h*dv
t+=0.5*h
plot(t, x_half,'r+')
#proceed into stepping loop
i=1
while i<nsteps:
   dx, dv = f(x_half, v_half)
   #update from t to t + h
   x_full += h*dx
   v_full += h*dv
   t+=0.5*h
   plot(t, x_full,'b+-')
   t+=0.5*h
   #update from t+h/2 to t+3h/2
   #the following is the full step tendency
   dx,dv = f(x_full,v_full)
   x_half += h*dx
   v_half += h*dv
   t+=0.5*h
   plot(t, x_half,'r+')   
   i+=1
show()