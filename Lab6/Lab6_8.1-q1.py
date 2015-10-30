from numpy import arange
import matplotlib.pyplot as plt

RC1 = 0.1
RC2 = 0.01
RC3 = 1.0
RC = RC3

def V_in(t):
	
	if int(2*t)%2==0:
		return 1
	else:
		return -1

def f(x,t):
    return (V_in(t) - x)/RC

a = 0.0
b = 10.0
N = 1000
h = (b-a)/N

tpoints = arange(a,b,h)
xpoints = []
x = 0.0

for t in tpoints:
    xpoints.append(x)
    k1 = h*f(x,t)
    k2 = h*f(x+0.5*k1,t+0.5*h)
    k3 = h*f(x+0.5*k2,t+0.5*h)
    k4 = h*f(x+k3,t+h)
    x += (k1+2*k2+2*k3+k4)/6

plt.plot(tpoints,xpoints)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Low pass filter for RC = " + str(RC))
plt.show()
