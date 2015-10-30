Eric Yeung
999979974

#IMPORT modules division, math, pylab, numpy, and pyplot

#DEFINE spring constant k

#DEFINE mass of the relativistic particle m

#DEFINE tstep, initial and final times

#CREATE and DEFINE time array

#Initialize arrays x and v for position and velocities

#OBTAIN intial conditions x_0 and v_0, make these the first entrys in arrays

#Intialize loop variable i=0

#While i< tlen-1:
	#Calculate x[i+1] and v[i+1] from x[i] and v[i] using euler method
 	   #x[i+1] = x[i] + v[i]*dt
	   #v[i+1] = v[i] - k/m*x[i]*(1- (v[i])**2/c**2)**(3/2)*dt
	#INCREMENT i
	#Loop repeats until i < tlen-1 is no longer satsified 

#Intialize new non-relativistic arrays xn and vn for NR positions and NR velocities

#OBTAIN intial conditions xn_0 and vn_0, make these the first entrys in arrays

#Intialize loop variable j=0, use a different index

#While j< tlen-1:
	#Calculate xn[j+1] and vn[j+1] from xn[j] and vn[j] using euler method
    	#xn[j+1] = xn[j] + vn[j]*dt

    	#REMOVE the relativistic term in the differential equation for velocity
    	#vn[j+1] = vn[j] - k/m*xn[j]*dt
	#INCREMENT j
	#Loop repeats until j < tlen-1 is no longer satsified 

#SHOW explicitly the difference after time evolution by comparing last entrys in arrays

#PLOT subplots comparing x vs t and xn vs t

#PLOT subplots comparing v vs t and vn vs t
