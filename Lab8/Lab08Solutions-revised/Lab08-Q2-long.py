#gravity_wave.py
#By Paul Kushner for PHY407
#Generate the displacement of the free surface under the action of gravity in an infinitely deep layer

#Solve Laplace equation \nabla^2\phi=0 by iterative relaxation subject to the boundary conditions that
#i) phi = 0 for all boundary
#ii) The top and bottom values of $\psi$ are fixed for all $x$.
#The bottom boundary condition is actually incorrect for the shallow water wave problem. In this simulation we impose the Dirichlet boundary condition $\phi=0$ when we should apply the neumann condition $\partial\phi/\partial z=0$ at the bottom $z=-H$.

#import necesssary routines
from numpy import copy, linspace, empty, pi, exp, zeros, meshgrid, arange
from pylab import plot, contourf, colorbar, xlim, ylim, xlabel, ylabel,clf, title,\
     pause, figure, show, contour, subplot, gradient, axis, quiver, legend, savefig
from time import time

from pickle import dump

#Define routine to calculate Laplacian
def laplace(phi,niter=1000):
   #We require a temporary array for this
   #create a copy
   phicopy = copy(phi)
   #define target
   target = 5e-5 #m^2/s
   #keep iterating until target is reached or until failure
   #for each iteration
   for i in range(niter):
      #Do nearest neighbour averaging in interior domain away from boundaries
      phi[1:-1,1:-1]=(phi[:-2,1:-1]+phi[2:,1:-1]+phi[1:-1,:-2]+phi[1:-1,2:])/4 
      #calculate difference between current and previous
      diff = abs(phi-phicopy)
      #if difference is within target then break out of this loop
      if (diff.max()<target):
         print 'iteration i',i
         print 'diff', diff.max()
         return i
         break
      phicopy= copy(phi)
      #if not then just continue
      
   #if convergence isn't reached after the number of iterations is exceeded then quit
   if (diff.max()>target):
      print 'not converged'
      exit()
         

#The tendencies are at the surface
#$\partial\phi(x,0,t)/\partial t=-g\eta$ and $\partial\eta(x,t)/\partial t=\partial\phi(x,0,t)/\partial z$.
#The latter expression requires a finite difference calculation implemented below
def dfdt(phi,eta,dz):
   #implement w = deta/dt=dphi/dz, dphi/dt = -g*eta
   #the phiz is a finite difference involving surface quantities.
   phiz=(phi[-1,:]-phi[-2,:])/dz
   return -g*eta,phiz

#Geometry, parameters
#define density, gravity
rho = 1.0e3
g=9.8

#define domain parameters
L=800.0
D=50.0
#define numerical parameters
h = 0.01 #seconds
#define initial time
t=0
#define nsteps 
nsteps = 4100
#number of snapshots
snap_int = 25
nsnaps = nsteps/snap_int
#define nsteps in x and z
npts_x   = 1601
npts_z   = 101

output = open("Lab08-Q2d-h=%f-nsteps=%i-nx=%i-nz=%i.p"%(h,nsteps,npts_x,npts_z),"wb")
#dump((L,D,h,nsteps,snap_int,nsnaps,npts_x,npts_z),output)


#define eta(x,t), phi(x,z), surface phi0(x,t)
#the eta_array, phi_array are defined for all time points
eta_array = empty([nsnaps,npts_x],float) #in m
phi_array = empty([nsnaps,npts_z,npts_x],float) #in m^2/s
#define time axis
time_array = empty([nsnaps],float) #in s
time_array[0] = 0.
#define grid of points in x and z
x=linspace(0,L,npts_x)
z=linspace(-D,0,npts_z)
#define resolution in x and z
dx = x[1]-x[0] #spatial grid resolution
dz = z[1]-z[0]

#width of pulse
Delta = 2.0#m
#height of pulse
A = -1.0
#initialize $\eta$ 
#find initial time
clockstart = time()
#initial eta (on full timestep)
#eta_full is the full timestep version of eta
eta_full = A*exp(-(x-L/2.0)**2/(Delta)**2)

#set up eta1, eta2 as zero arrays.
#eta_half is 
eta_half = empty(npts_x, float)

#phi0 is the surface potential, which is updated
#Set initial phi to zero indicating displaced interface and state of rest.
phi0_full = zeros(npts_x, float)
phi0_half = zeros(npts_x, float)

#Initialize phi to zero. The boundary conditions are zero for all times. 
phi = zeros((npts_z, npts_x), float)



#initial time step 
#variables that will be time stepped are eta, phi0
#need separate eta_0 and phi_0
#simple euler forward for first half step
dphi0, deta =dfdt(phi,eta_full,dz)
eta_half  = eta_full+0.5*h*deta
phi0_half = phi0_full+0.5*h*dphi0

#phi01 is surface condition
phi[-1,:]=phi0_half #update the surface condition

#solve for laplacian(phi)=0
#solve for #phi field consistent with this
laplace(phi)

#proceed into stepping loop
#for each time step
iterpoints = []
snap_num = 0 #snapshot number
eta_array[snap_num,:]   = eta_full[:]
phi_array[snap_num,:,:] = phi
time_array[snap_num] = t
snap_num += 1


for i in range(1,nsteps):
   loop_time = time()
   #update from t to t + h
   dphi0, deta =dfdt(phi,eta_half,dz)
   #update from t to t + h
   eta_full[:] += h*deta
   phi0_full[:]+= h*dphi0
   #update surface condition at t+h
   phi[-1,:] = phi0_full
   #solve phi consistent with this
   #and save number of iterations for plotting
   iterpoints.append(laplace(phi))
   #update time
   t+=h
   #save diagnostics
   if i%snap_int==0:
      eta_array[snap_num,:]   = eta_full[:]
      phi_array[snap_num,:,:] = phi
      time_array[snap_num] = t
      snap_num += 1
   #update from t+h/2 to t+3h/2
   #the following is the full step dfdt
   dphi0,deta = dfdt(phi,eta_full,dz)
   eta_half[:] += h*deta
   phi0_half[:] += h*dphi0
   phi[-1,:] = phi0_half[:]
   laplace(phi)
   #plot
   clf()
   plot(x,eta_full,'k',x,eta_half,'r')
   ylim([-2*A,2*A])
   xlim([L/2-L/5,L/2+L/5])
   title('t=%3.3f'%t)
   pause(0.01)
   i+=1
   print 'time in this loop',  time() - loop_time
total_time = time()-clockstart
print 'time:' , total_time

dump((eta_array, phi_array, x, z,time_array),output)
output.close()