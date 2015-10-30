# Script demonstrating some examples of formatting floats and integers into strings
# Oliver Watt-Meyer
# October 21, 2015

from numpy import pi

a=pi

# formatting a float (round to 2nd decimal place)
print a
print 'a = %.2f' % a

# formatting with scientific notation
b=1.836e-11
print 'b = %.2e' % b

# printing an integer
c=8
print 'c = %d' % c

# formatting multiple numbers in one string
print 'a = %.2f, b=%.1e, c=%d' % (a,b,c)

# padding an integer with zeros
for i in range(5):
    print 'i = %03d' % i
    
# padding an integer with spaces
for i in range(5):
    print 'i = %3d' % i   