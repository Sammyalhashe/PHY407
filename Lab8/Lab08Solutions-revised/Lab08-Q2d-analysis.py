#Lab08-Q2d-analysis
#By Paul Kushner for PHY407
#Analyze solutions from previous integration of gravity wave code


#import necesssary routines
from numpy import copy, linspace, empty, pi, exp, zeros, meshgrid, arange, concatenate
from pylab import plot, contourf, colorbar, xlim, ylim, xlabel, ylabel,clf, title,\
     pause, figure, show, contour, subplot, gradient, axis, quiver, legend, savefig,grid, text
from time import time

import pickle
input = open("Lab08-Q2d-h=0.010000-nsteps=4100-nx=1601-nz=101.p","rb")
(eta_array, phi_array, x, z,time_array) = pickle.load(input)
levs= concatenate([arange(-1,0.0,0.05),arange(0.05,1.05,0.05)])
figure(1,figsize=(10,8))
clf()
contourf(x,time_array,eta_array,levs)
contour(x,time_array,eta_array,levs,colors='k')
L=x[-1]
xlim((L/2-L/10,L/2+L/10))
xlabel('x(m)')
ylabel('t (s)')
show()
savefig('Lab08-Q2d-Figure1.pdf')
D=-z[0]
times_index = [0,8,40,160]
times_string = ['0','2','10','40']
figure(2,figsize=(12,10))
clf()
subplot_order = [1,3,5,7,2,4,6,8]
for j in range(len(times_index)):
   subplot(4,2,subplot_order[2*j])
   plot(x,eta_array[times_index[j],:])
   xlim((L/2-L/5,L/2+L/5))
   place_text_y = eta_array[times_index[j],:].max()*0.8+eta_array[times_index[j],:].min()*0.2
   place_text_x = 250
   text(place_text_x,place_text_y,'t=%3.2f s'%time_array[times_index[j]])
   grid('on')
   subplot(4,2,subplot_order[2*j+1])
   contour(x,z,phi_array[times_index[j],:,:],colors='k')  
   contourf(x,z,phi_array[times_index[j],:,:])
   colorbar(orientation='horizontal')
   xlim((L/2-L/5,L/2+L/5)) 
   ylim((-D/10,0))
   grid('on')
show()
savefig('Lab08-Q2d-Figure2.pdf')
figure(3)
clf()
h=0.25
#plot profile at a few different times
times_index = [int(30/h),int(31/h),int(32/h),int(33/h)]
for j in range(len(times_index)):
   plot(x,eta_array[times_index[j],:]+0.5*j)
   
xlim((L/2,L/2+L/5)) 
grid('on')
legend(('%3.2f'%time_array[times_index[0]],'%3.2f'%time_array[times_index[1]],'%3.2f'%time_array[times_index[2]],'%3.2f'%time_array[times_index[3]]),loc='lower right')
show()

#the following was taken from the figure above:
g=9.81
#picked up out some features in the wave packet to track phase velocity
t1,x11,x12,x13=(30.0,454.0,469.0,491.0)
t2,x21,x22,x23=(33.0,467.0,485.0,511.0)
xpos = [x11,x12,x13,x21,x22,x23]
ypos = [0.1,0.1,0.1,1.6,1.6,1.6]
strs = ['k.','kx','k+','k.','kx','k+']
for j in range(6):
   plot(xpos[j],ypos[j],strs[j])

delta_t = t2-t1
delta_x1 = x21-x11
delta_x2 = x22-x12
delta_x3 = x23-x13
lam_11,lam_12 = (x12-x11,x13-x12)
lam_21,lam_22 = (x22-x21,x23-x22)
lam_1_mean = 0.5*(lam_11+lam_21)
lam_2_mean = 0.5*(lam_12+lam_22)
c_1 = delta_x1/delta_t
c_2 = delta_x2/delta_t
c_3 = delta_x3/delta_t
c_1_est = 0.5*(c_1+c_2)
c_2_est = 0.5*(c_2+c_3)
c_disp_1 = (g*lam_1_mean/(2*pi))**0.5
c_disp_2 = (g*lam_2_mean/(2*pi))**0.5
pdiff_1 = abs(c_1_est-c_disp_1)/c_disp_1*100
pdiff_2 = abs(c_2_est-c_disp_2)/c_disp_2*100

#group velocity calculation
#picked out peak amplitudes in wave packet for group velocity
t1,x11g,x12g,x13g=(30.0,443.,448.,454.)
t2,x21g,x22g,x23g=(33.0,449.,454.,460.)
x1g = (x11g+x12g+x13g)/3.0
x2g = (x21g+x22g+x23g)/3.0
cg_est = (x2g-x1g)/(t2-t1)
lam_group_1 = x13g-x11g
lam_group_2 = x23g-x21g
lam_group_est = 0.5*(lam_group_1+lam_group_2)
cg_from_disp = (g*lam_group_est/(2*pi))**0.5/2
pdiff_group_calc = abs(cg_est-cg_from_disp)/cg_est*100
xpos = [x11g,x12g,x13g,x21g,x22g,x23g]
ypos = [0,0,0,1.5,1.5,1.5]
strs = ['r.','rx','r+','r.','rx','r+']
for j in range(6):
   plot(xpos[j],ypos[j],strs[j])
show()
print 'c_1_est %3.2f m/s, c_disp_1 %3.2f m/s, percent diff %3.0f'%(c_1_est, c_disp_1, pdiff_1)
print 'c_2_est %3.2f m/s, c_disp_2 %3.2f m/s, percent diff %3.0f'%(c_2_est, c_disp_2, pdiff_2) 

print 'estimated group velocity, group velocity estimated from dispersion relation, pc diff: %3.2f %3.2f %3.0f'%(cg_est,cg_from_disp,pdiff_group_calc)
savefig('Lab08-Q2d-Fig3.pdf')
xlim((L/2,L/2+L/20.0))
ylim((-0.1,0.1))
show()
savefig('Lab08-Q2d-Fig3-zoom.pdf')
