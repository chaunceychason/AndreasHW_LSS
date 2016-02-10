from sys import exit
import os
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt

#import numpy as np
#import pandas as pd
#from scipy import stats  
#from matplotlib import rc
#import matplotlib.pyplot as plt


def plot_basic( xlist, ylist, title, xlab, ylab, psize, yflip):
	print("Entered Basic Plot Function")
	plot_title=" Blank Title "   
	x_axis="Blank X"
	y_axis="Blank Y"
	pointsize = 5
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#Choose which type of plot you would like: Commented out.
	
	
	if True:
		plot_title = title
		x_axis = xlab
		y_axis = ylab
 		pointsize = psize

	
	#plt.plot(xlist, ylist)
	#plt.plot(xlist, ylist, linestyle="", marker="o", markeredgecolor=None, markeredgewidth=0.0)
	#plt.plot(scale1, Vmaxarray, linestyle="-", marker="o")
	plt.scatter(xlist, ylist, s=pointsize, lw=0)
	
	plt.title(plot_title)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	#plt.yscale("log")
	
	if yflip == True:
		plt.ylim(max(ylist), min(ylist))

	#plt.savefig(figure_name)
	#print("Saving plot: %s" % figure_name)
	print
	
	#Comment out to over plot curves.			
	#plt.clf()

	"""
	clearmass  = []
	clearscale = []
	clearVmax  = []
	return clearmass, clearscale, clearVmax
	"""
	dummy = 1
	return dummy

def checkcondition(condition1, condition2):
	if (condition1 or condition2) == False:
		print("Conditions evaluated to FALSE! This should break. ")
		return False
	else:
		#The conditions are true!
		return True

"""
================================================================
================================================================

			    7MM766Yb.  7MM766Mq.           
			    MM    `Yb. MM   `MM.          
			    MM     `Mb MM   ,M9 M******A' 
			    MM      MM MMmmdM9  Y     A'  
			    MM     ,MP MM  YM.       A'   
			    MM    ,dP' MM   `Mb.    A'    
			  .JMMmmmdP' .JMML. .JMM.  A'     
			                          A'      
			                         A'       
HOMEWORK #2.  Written by:    Nicholas Chason 
----------------------------------------------------------------
INSTRUCTIONS: 
A.  - i.  Make a plot of DEC vs. RA showing positions. 
	- ii. Make a plot of Mr vs z showing all galaxies. (invert y-ax)

B.  - i.  Make a plot of g-r color distribtuion of galaxies. 
	- ii. Using a cut between red and blue of g-r=0.75 calculate the 
			blue fraction in the sample.
C.  - i.  Plot the r-band luminosity function of galaxies: 
			log(dn/dMr) vs. Mr. [units: dn {h-1 Mpc}] in 0.1 M bins. 
			Use a median z survey depth of 0.1. 
			DR7 coverage:  7675.2 deg^2.  (2.295 steradians)

D.  - i.  Construct three volume limited sub-samples of the full data-
			set containing galaxies more luminous than -20, -19, & -18. 
			List the z-bounds, volume, and number of galaxies in each. 
	- ii. Calculate the new blue galaxy fractions for these samples. 

	- iii. Go to E. for plot instructions. 
E.  - i.  Make a single plot of the luminosity function meausured from
			the three volume-limited samples from D. 

F.  - i. Make a plot of the luminosity function using the full flux-
			limited sample, applying a 1/V_max weight. How does the
			weighting change the LF? 
================================================================
================================================================
"""
#Sets the datafile name for opening
datafilename = "./SDSS_DR7.dat"

"""
DATA FILE STRUCTURE: SDSS_DR7.dat 
----------------------------------------------------------------
COLUMNS {0, 1, 2, 3, 4} x ROWS {550166} : 

RA, DEC, z, Mg, Mr 

Note: Each row contains infomation for a single galaxy. 	
----------------------------------------------------------------
"""

if os.path.exists(datafilename):
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")    

condition1 = True
condition2 = True

with open(datafilename) as fp:

	#Initializes the arrays for plotting
	RA_LIST   = []      
	DEC_LIST  = []		
	z_LIST    = []		
	abs_g_mag_LIST  = []	
	abs_r_mag_LIST  = []	

	for line in fp:
		#print line
		if checkcondition(condition1, condition2) == False:
			print("Condition FALSE.")
			break
		else:
			pass		#Do nothing. Continue on with code. 


		splitline = line.split()  #Divide the non-hashed data into splitline.

		RA_value   = float(splitline[0])        #Sets RA
		DEC_value  = float(splitline[1])		#Sets Dec
		z_value    = float(splitline[2])		#Sets  z
		abs_g_mag  = float(splitline[3])		#Sets  g-band mag
		abs_r_mag  = float(splitline[4])		#Sets  r-band mag

		RA_LIST.append(RA_value)
		DEC_LIST.append(DEC_value)
		z_LIST.append(z_value)
		abs_g_mag_LIST.append(abs_g_mag)
		abs_r_mag_LIST.append(abs_r_mag)


#PLOT FUNCTION. { X,  Y}
# A.
#======================================================
title_label = "title string"
x_label = "x axis"
y_label = "y axis"
x_data  = RA_LIST
y_data  = DEC_LIST 
pointsize = 1
yflip = False

if len(x_data) != len(y_data):
	print("ERROR! X and Y DATA LENGTHS ARE DIFFERENT!")
	print("Length: x_data: %g" % len(x_data))
	print("Length: y_data: %g" % len(y_data))
else:
	print("Length: x_data: %g" % len(x_data))
	print("Length: y_data: %g" % len(y_data))

dummyvariable = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip)
figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()
#======================================================
# A. part II
#======================================================
title_label = "title string"
x_label = "x axisAII"
y_label = "y axisAII"
x_data  = z_LIST
y_data  = abs_r_mag_LIST
pointsize = 1
yflip = True

if len(x_data) != len(y_data):
	print("ERROR! X and Y DATA LENGTHS ARE DIFFERENT!")
	print("Length: x_data: %g" % len(x_data))
	print("Length: y_data: %g" % len(y_data))
else:
	print("Length: x_data: %g" % len(x_data))
	print("Length: y_data: %g" % len(y_data))

dummyvariable = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip)
figure_name=os.path.expanduser('~/Feb9LSSHW2_plot2_A' +'.png')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
#======================================================


