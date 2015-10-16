#PHY407, Fall 2015: Lab01, Q2
#Author: DUONG, BANG CHI

#Relativistic particle on a spring.

#Use Euler Cromer method to do the integration (more stable than Euler).

#---------Import statements-----------
#import zeros from numpy
from numpy import zeros, sqrt, linspace, argmin
from numpy.fft import rfft
from pylab import plot, show, figure, subplot, title, xlabel, ylabel, grid, legend, suptitle, savefig, xlim, axvline, tight_layout
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
T = 120.0 #s

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
ft_vRel_1m = rfft(vRel_1m)
#xc case
x0 = xc #m
v0 = 0.0 #m/s
xClass_xc, vClass_xc, dummy = RHS_Class(k, m, x0, v0, dt, T)
xRel_xc, vRel_xc, dummy = RHS_Relativistic(k, m, c, x0, v0, dt, T)
ft_xc = rfft(xRel_xc)
ft_vRel_xc = rfft(vRel_xc)
#10xc case
x0 = 10*xc #m
v0 = 0.0 #m/s
xClass_10xc, vClass_10xc, dummy = RHS_Class(k, m, x0, v0, dt, T)
xRel_10xc, vRel_10xc, dummy = RHS_Relativistic(k, m, c, x0, v0, dt, T)
ft_10xc = rfft(xRel_10xc)
ft_vRel_10xc = rfft(vRel_10xc)
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
    tight_layout()
    savefig('Lab03-Q1a-Fig'+str(j+1)+'.png')
#------------------------------------------------------------------------
# Relativistic period
x0_range = linspace(1,10*xc,1000)

def g(x):
    return c*sqrt(((0.5*k*(x0**2-x**2)*(2.0*m*c**2+0.5*k*(x0**2-x**2))) \
                      /(m*c**2+0.5*k*(x0**2-x**2))**2))

def integrand(x):
    return 4.0*g(x)**-1

T_relativistic = zeros(len(x0_range))
# x and w in interval -1 to 1
x,w = gaussxw(200)
for i in range(len(x0_range)):
    s = 0.0
    x0 = x0_range[i]
    xp = 0.5*x0*x + 0.5*x0
    wp = 0.5*x0*w
    
    for j in range(200):
        s += wp[j]*integrand(xp[j])
    T_relativistic[i] = s

# Plots for Fourier Transformation, and the method used in lab3-Q2b
figure()
subplot(211)
plot(abs(ft_1m)/max(abs(ft_1m)),label='pos_1m')
plot(abs(ft_xc)/max(abs(ft_xc)),label='pos_xc')
plot(abs(ft_10xc)/max(abs(ft_10xc)),label='pos_10xc')
axvline(120/T_relativistic[0], label='freq_1m', color = 'b', linestyle = '--')
axvline(120/T_relativistic[argmin(abs(x0_range-xc))], label='freq_xc', color = 'g' , linestyle = '--')
axvline(120/T_relativistic[999], label='freq_10xc', color = 'r', linestyle = '--')
xlabel('Frequency f')
ylabel('$|\^x(f)|$ / $|\^x(f)|_{max}$')
legend(loc='upper right')
grid('on')
xlim(0,100)
tight_layout()
savefig('Lab03-Q1a-Position.png')

subplot(212)
plot(abs(ft_vRel_1m)/max(abs(ft_vRel_1m)),label='vel_1m')
plot(abs(ft_vRel_xc)/max(abs(ft_vRel_xc)),label='vel_xc')
plot(abs(ft_vRel_10xc)/max(abs(ft_vRel_10xc)),label='vel_10xc')
axvline(120/T_relativistic[0], label='freq_1m', color = 'b', linestyle = '--')
axvline(120/T_relativistic[argmin(abs(x0_range-xc))], label='freq_xc', color = 'g' , linestyle = '--')
axvline(120/T_relativistic[999], label='freq_10xc', color = 'r', linestyle = '--')
xlabel('Frequency f')
ylabel('$|\^v(f)|$ / $|\^v(f)|_{max}$')
legend(loc='upper right')
grid('on')
xlim(0,100)
tight_layout()
savefig('Lab03-Q1a-Velocity.png')

show()

# Characteristic Frequencys
print 120/T_relativistic[0], 120/T_relativistic[argmin(abs(x0_range-xc))], 120/T_relativistic[999]

# Periods
print T_relativistic[0], T_relativistic[argmin(abs(x0_range-xc))], T_relativistic[999]
