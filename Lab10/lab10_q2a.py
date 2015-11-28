# PHY407, Fall 2015, Question 2a
# Author: DUONG, BANG CHI

# Import modules
from __future__ import division, print_function
from pylab import plot, imshow, gray, draw, savefig, title, \
     xlabel, ylabel, figure, xlim, ylim, axvline, axhline, show,\
     ion, zeros, arange

from visual import sphere, rate, display
from random import randint

# Define constants
L = 101 # length  of the lattice

i = 50
j = 50

i_list = []
j_list = []

# Animation

d = display()
s = sphere()
s.pos = i,j,0

d.autoscale = False
for t in arange(1e4):
	#rate(30)
	s.pos = i-50,j-50,0
	
	direction = randint(1,4)	
	if direction == 1: 
		if i==L: continue
		i+=1 # move right
	elif direction == 2:
		if i==0: continue
		i-=1 # move left
	elif direction == 3:
		if j==L: continue
		j+=1 # move up
	elif direction == 4:
		if j==0: continue
		j-=1 # move down

	i_list.append(i)
	j_list.append(j)

# Plot the path of the particle
plot(i_list,j_list)
title('Path travelled by the particle')
xlabel('x position')
ylabel('y position')
xlim(-10,L+10)
ylim(-10,L+10)
axvline(-1); axvline(L)
axhline(-1); axhline(L)
show()
		


