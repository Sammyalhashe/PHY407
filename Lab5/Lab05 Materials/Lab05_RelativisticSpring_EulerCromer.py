#PHY407, Fall 2015: Lab01, Q2
#Author: Paul Kushner
#Relativistic particle on a spring.
#Integrate dp/dt = F = -kx using Euler forward time stepping.
#Carry out and compare classical and relativistic calculations.
#Use Euler Cromer method to do the integration (more stable than Euler).

#---------Import statements-----------
#import zeros from numpy
from numpy import zeros, pi, dot, arange
from numpy.fft import rfft
from pylab import plot, show, figure, subplot, title, xlabel, legend, suptitle, savefig, xlim
from gaussxw import gaussxw



#weights and sampling points for each case

#define function to calculate rhs in dx/dt = f for classical calculation
#input: initial conditions, length of integration, time step, parameters k, m, c
#output: solution position and velocity
def RHS_Class(k, m, x0, v0, dt, T):
    #number of steps is ratio of total time to timestep
    N_Steps = int(T/dt)
    #initialize arrays: u1, u2, t as zero
    u1 = zeros(N_Steps)
    u2 = zeros(N_Steps)
    t  = zeros(N_Steps)
    #initialize u1[0] as input x0
    #initialize u2[0] as input v0    
    u1[0] = x0
    u2[0] = v0
    t[0] = 0
    #for each time step
    
    for j in range(N_Steps-1):
        #R_1 = u2
        #R_2 = -k*u1/m

        r2 = -k*u1[j]/m 
        u2[j+1] = u2[j] + r2*dt        
        #use Euler-Cromer   
        r1 = u2[j+1]
        u1[j+1] = u1[j] + r1*dt
        t[j+1] = dt*(j+1)
    return u1, u2, t

#define function to calculate rhs in dx/dt = f for relativistic calculation
#input: initial conditions, length of integration, time step, parameters k, m, c
#output: solution position and velocity
def RHS_Relativistic(k, m, c, x0, v0, dt, T):
    #number of steps is ratio of total time to timestep
    N_Steps = int(T/dt)
    #initialize arrays: u1, u2, t as zero
    u1 = zeros(N_Steps)
    u2 = zeros(N_Steps)
    t  = zeros(N_Steps)
    #initialize u1[0] as input x0
    #initialize u2[0] as input v0    
    u1[0] = x0
    u2[0] = v0
    t[0] = 0
    #for each time step
    
    for j in range(N_Steps-1):
        r2 = -k*u1[j]*(1-u2[j]**2/c**2)/m   
        u2[j+1] = u2[j] + r2*dt
        #Euler-Cromer
        r1 = u2[j+1]
        u1[j+1] = u1[j] + r1*dt
        t[j+1] = dt*(j+1)
    return u1, u2, t
  
#---------Main Program-----------

#set parameters
k = 12.0 #N/m
m = 1.0 #kg
c = 3.0e8 #m/s

#time step and time elapsed
dt = 1e-4 #s
T = 10.0 #s

#determine xc as per lab question.
xc = c*(m/k)**0.5
print 'xc is', xc/1e6, 'thousands of kilometers'


#call calculation and return x, v, t for each case
#1m case
x0 = 1.0 #m
v0 = 0.0 #m/s
xClass_1m, vClass_1m, t = RHS_Class(k, m, x0, v0, dt, T)
xRel_1m, vRel_1m, dummy = RHS_Relativistic(k, m, c, x0, v0, dt, T)
ft_1m = rfft(xRel_1m)
#xc case
x0 = xc #m
v0 = 0.0 #m/s
xClass_xc, vClass_xc, dummy = RHS_Class(k, m, x0, v0, dt, T)#
xRel_xc, vRel_xc, dummy = RHS_Relativistic(k, m, c, x0, v0, dt, T)
ft_xc = rfft(xRel_xc)
#10xc case
x0 = 10*xc #m
v0 = 0.0 #m/s
xClass_10xc, vClass_10xc, dummy = RHS_Class(k, m, x0, v0, dt, T)
xRel_10xc, vRel_10xc, dummy = RHS_Relativistic(k, m, c, x0, v0, dt, T)
ft_10xc = rfft(xRel_10xc)
#Plot each of the cases

for j in range(3):
    #choose the case
    if j==0:
        xClass, xRel, vClass, vRel = (xClass_1m, xRel_1m, vClass_1m, vRel_1m)
        x0 = 1.0
    if j==1:
        xClass, xRel, vClass, vRel = (xClass_xc, xRel_xc, vClass_xc, vRel_xc)
        x0 = xc
    if j==2:
        xClass, xRel, vClass, vRel = (xClass_10xc, xRel_10xc, vClass_10xc, vRel_10xc)
        x0 = 10*xc
    
        
    #plot the two calculations as an overlay, for both x and v
    figure(j)
    
    subplot(2,1,1)
    plot(t, xRel, t, xClass)
    legend(('Relativistic', 'Classical'))
    xlabel('time(s)')

    subplot(2,1,2)
    plot(t, vRel, t, vClass)
    xlabel('time(s)')
    title('v')
    suptitle('x0 = '+str(x0)+'m')
    savefig('Lab01-Q2-Fig'+str(j+1)+'.pdf')


show()
