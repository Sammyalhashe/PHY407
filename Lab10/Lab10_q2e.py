#################################################################
# This program simulates diffusion limited aggregation on an LxL grid.
# Particles are initiated until the centre point is filled.
#################################################################
from random import randrange, randint
from numpy import arange,zeros,log,sqrt,nonzero
from scipy.optimize import curve_fit
from matplotlib.pyplot import plot,scatter,draw,savefig,title,legend,\
     xlabel,ylabel,figure,xlim,ylim,ion,xscale,yscale,show

def nextmove(x,y):
    direction=randrange(1,5)

    # randomly choose a direction
    # 1 = up, 2 = down, 3 = left, 4=right

    if direction==1:
        #move up
        y+=1
    elif direction==2:
        #move down
        y-=1
    elif direction==3:
        #move right
        x+=1 
    elif direction==4:
        #move left
        x-=1
    else:
        print "error: direction isn't 1-4"
        
    return x,y

#############################################################
# main program starts here
#############################################################

#ion()

Lp=201 # size of domain
N=2000 # number of particles
anchored=zeros((Lp,Lp),dtype=int) # array to represent whether each gridpoint has an anchored particle
anchored_points=[[],[]] # list to represent x and y positions of anchored points
centre_point=(Lp-1)/2 # middle point of domain

# set up animation of anchored points
#animation_interval=50 # how many moves to make before updating plot of Brownian motion
##figure(1)
##title('Diffusion-limited aggregation run for 2000 particles') 
###moving_plot=plot(centre_point,centre_point,'.r',markersize=10)
##anchored_plot=plot(anchored_points[0],anchored_points[1],'.k',markersize=10)
##xlim([-1,Lp])
##ylim([-1,Lp])
##xlabel('x')
##ylabel('y')

for j in range(N):
    xp = randint(0,Lp)
    yp = randint(0,Lp)
    i=0 # counter to keep track of animation of moving particle
    
    not_stuck=True
    while not_stuck:

        if xp == centre_point and yp == centre_point: # Check if the particle has reached the centre
            anchored[xp,yp]=1
            anchored_points[0].append(xp)
            anchored_points[1].append(yp)
##            anchored_plot[0].set_xdata(anchored_points[0])
##            anchored_plot[0].set_ydata(anchored_points[1])
##            draw()

            not_stuck=False
            
        elif anchored[xp-1:xp+2,yp-1:yp+2].any(): # Check if particle is adjacent to an anchored particle
            anchored[xp,yp]=1
            anchored_points[0].append(xp)
            anchored_points[1].append(yp)
##            anchored_plot[0].set_xdata(anchored_points[0])
##            anchored_plot[0].set_ydata(anchored_points[1])
##            draw()
            
            not_stuck=False
        
        elif xp==0 or xp==Lp-1 or yp==0 or yp==Lp-1:
            if xp==0:
                xp==Lp-1
            elif xp==Lp-1:
                xp==0
            elif yp==0:
                yp==Lp-1
            elif yp==Lp-1:
                yp==0
##            draw()

            not_stuck=False

        else: # If neither of the above, move particle and continue while loop
            i+=1
            xp,yp=nextmove(xp,yp)

        #print i, xp, yp

#savefig('lab10_q2d.jpg')

###------------------------Part e-----------------------------###

# Calculate fractal dimension
    # Define arrays for length and count
length = arange(2,Lp+1)
count = zeros(len(length))

    # Count the 'no. of anchored points' that fall inside the box of length l
tracking = 0
for i in range(len(anchored_points[0])):
    tracking += 1
    for j in range(len(length)):
        if (Lp+1-length[j])/2 < anchored_points[0][i] < (Lp+1+length[j])/2 \
           and (Lp+1-length[j])/2 < anchored_points[1][i] < (Lp+1+length[j])/2:
            count[j] += 1
    print tracking

    # Do a log transformation on length and mass
length_log = log(length)
count_log = log(count)

length_log_linear = length_log[nonzero(length_log < 4.)]
count_log_linear = count_log[nonzero(length_log < 4.)]

    # Do Linear Fit
def LinFit(x_value,*p):
    return p[0] + p[1] * x_value

FitPara, FitCov = curve_fit(LinFit, length_log_linear, count_log_linear, 
                            p0=(0,1), maxfev = 100000)
Fit_count_log_linear = LinFit(length_log_linear,*FitPara)

    # Analysis
print 'With Lp = 201, and N = 2000 particles, the estimate Fractal dimension k = {} +/- {}'.format(FitPara[1],sqrt(FitCov[1,1]))

    # Plots
figure()
plot(length, count)
xscale('log')
yscale('log')
xlabel('Length of the box (log scale)')
ylabel('Number of anchored points (log scale)')

figure()
scatter(length_log, count_log)
plot(length_log_linear, Fit_count_log_linear, 'r', label = 'Linear fit')
xlabel('Log(Length of the box)')
ylabel('Log(Number of anchored points)')
legend().draggable()

show()
