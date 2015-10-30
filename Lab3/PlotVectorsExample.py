#plot a polynomial function of x and y and its gradient
#import required functions
from pylab import quiver, figure, contour, axis, xlabel, ylabel

from numpy import sin,ones, pi, linspace, zeros

#adjust vec_spacing until you get the vectors to be placed at the right interval to your liking
vec_spacing = 1 #Try 4 or 5
#x domain
nx = 51
x = linspace(-1,1,nx)
#ydomain
ny = 101
y = linspace(-2,2,ny)

#length
f = zeros((ny,nx)) #flip x and y for plotting purposes
dfdx = zeros((ny, nx))
dfdy = zeros((ny, nx))

#The function is f = x**2+y**3, and its gradient is grad f = (2x, 3y^2)
for i in range(len(x)):
    for j in range(len(y)):
        f[j,i] = x[i]**2 + y[j]**3
        dfdx[j,i] = 2*x[i]
        dfdy[j,i] = 3*y[j]**2

figure(1)
#plot contours of f, using levels every 1.
contour(x,y,f,levels = linspace(-9,9,19),colors='k')

xlabel('x')
ylabel('y')

#plot vectors (gradient)
quiver(x[::vec_spacing],y[::vec_spacing],dfdx[::vec_spacing,::vec_spacing],dfdy[::vec_spacing,::vec_spacing])
#set up figure so that vectors are properly oriented.
axis('image')