from sys import exit
import os
import numpy as np
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt

#import numpy as np
#import pandas as pd
#from scipy import stats  
#from matplotlib import rc
#import matplotlib.pyplot as plt


def plot_basic( xlist, ylist, title, xlab, ylab, psize, yflip, pcounter):
	print("Entered Basic Plot Function")
	
	if len(xlist) != len(ylist):
		print("ERROR! X and Y DATA LENGTHS ARE DIFFERENT!")
		print("Length: x_data: %g" % len(x_data))
		print("Length: y_data: %g" % len(y_data))
	if len(xlist)==0 or len(ylist)==0:
		print("ERROR: list length is ZERO!")
		print("Length: x_data: %g" % len(x_data))
		print("Length: y_data: %g" % len(y_data))
	else:
		print("Length: x_data: %g" % len(x_data))
		print("Length: y_data: %g" % len(y_data))


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
		try:
			plt.ylim(max(ylist), min(ylist))
		except:
			print("uh.oh.... try except statement. check ylim.")
	#plt.savefig(figure_name)
	#print("Saving plot: %s" % figure_name)
	print
	
	figure_name=os.path.expanduser('~/Feb9LSSHW2_plot%s.png' % pcount)
	plt.savefig(figure_name)
	print("Saving plot: %s" % figure_name)
	plt.clf()

	dummy = pcounter + 1
	return dummy

def get_volume(max_z):
	#This volume uses max_redshift (max_z) to compute volume: 
	#	Note: Uses max_z, instead of median since it is a volume limited, instead of flux-lim. 
	"""
	#Note: 
		Could compute the volume accurately using the following resources:
  				http://home.fnal.gov/~gnedin/cc/    
		This site gives the distance between two redshifts. 

	#Alternative: Use the approximation that:  D ~= 3000 * z
	 This expression is actually more accurately D ~= 2950 * z. 	
	"""
	radius = 2950.*max_z 
	DR7_SAMPLE_COVERAGE_STERADIANS = 2.295   #Should be set to 2.295 sterad. (7675.2 deg^2)
	fract_vol = DR7_SAMPLE_COVERAGE_STERADIANS / (4.*np.pi)    #Fractional sky coverage DR7 to 4pi. 
	vol = fract_vol*(4./3.) * np.pi * (radius**3.) 
	return vol  #Units: [h-1 Mpc^3]

def get_bluefr(red, blue):
	if blue == 0:
		return 0
	blue_fraction =  blue / (1.0*red + 1.0*blue)
	return blue_fraction

def print_subsample_info(mag, x_data, red, blue ):
	print("==================\nM_r < -%d:  \n=====================" % abs(mag))
	print("The Redshift bounds   :[ %.4f --> %.4f ] " % (min(x_data) , max(x_data)))
	print("The Volume is         : %.4f [(h^-1 Mpc)^3]." % get_volume(max(x_data)))
	print("The # of Red Galaxies : %d " % red)
	print("The # of Blue Galaxies: %d " % blue)
	print("The Blue Fraction     : %.3f" % get_bluefr(red, blue))
	return True

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
	SDSS -	    MM      MM MMmmdM9  Y     A'  
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


NOTES:
	**********     NOTE: FIXED!   ***********
	# 
	#  USE:  {sort -k 4,4n filename > newfilename.} instead. 
	#
	#  TO FIND z:  sort -k4,4 -k3,3 SDSS_DR7.dat |  awk '{ if ($4 == "-19.00") print $0 }' 
	#
	#
	# * For part D, while trying to create a volume limited sample I 
	# 	attempted to sort the datafile by the 4th column (M_r) to 
	# 	speed up the data processing time. I did this using a bash 
	# 	command: { sort -k 4 SDSS_DR7.dat > SDSS_DR7sorted.dat } as
	# 	well as {sort -k 4,4 SDSS_DR7.dat > SDSS_DR7sorted.dat }. 
	# 	Both times, I noticed that it properly sorted for all but the
	# 	final two rows. They have been removed from the sample...
		
	# 	Lines Removed:
	# 	---------------------------------------------------- 
	# 	>	***	192.804726  47.070860 0.00183   -9.18 -11.41 
	# 	>	***	186.474874  33.539617 0.00104   -9.76 -10.00 

	* Next Issue... 	
================================================================
================================================================
"""
#-----------------------------------
#Sets the datafile name for opening
#-----------------------------------
#datafilename = "./SDSS_DR7.dat"
#datafilename = "./SDSS_DR7stable.dat"  #Sorted based on M_r (most negative last) 
#datafilename = "./SDSS_DR7condensed.dat"  #Sorted based on M_r (most negative last) 
datafilename = "./SDSS_DR7.dat"
ordered_file = True

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

#A counter to record the number of blue & red galaxies. g-r=0.75
g_r_more7p5counter, g_r_more7p5counter18, g_r_more7p5counter19, g_r_more7p5counter20 = 0, 0, 0, 0
g_r_less7p5counter, g_r_less7p5counter18, g_r_less7p5counter19, g_r_less7p5counter20 = 0, 0, 0, 0


#plot number counter
pcount = 1

#Variable that stores the index of the subsample cutoff point
current_index = 0
index_18subsample = 0
index_19subsample = 0
index_20subsample = 0

"""
=================================================
Sets the volume limiting redshift, from: 
=================================================
       RA          DEC        Z       M_r    M_g
-20: 203.645876  20.257367 0.16635  -20.00 -21.34
-19: 123.533843  48.703091 0.13640  -19.00 -20.95
-18: 57.696922   -5.360849 0.07342  -18.00 -19.16
"""
#z_lim20=0.16635 
z_lim20=0.15676
z_lim19=0.13640
z_lim18=0.07342

with open(datafilename) as fp:

	#Initializes the arrays for plotting
	RA_LIST   = []      
	DEC_LIST  = []		
	z_LIST    = []		
	abs_g_mag_LIST  = []	
	abs_r_mag_LIST  = []	
	gr_color_LIST   = []
	#n20subsample_indexes = []
	#n19subsample_indexes = []
	#n18subsample_indexes = []

	RA_LIST20   = []      
	DEC_LIST20  = []		
	z_LIST20    = []		
	abs_g_mag_LIST20  = []	
	abs_r_mag_LIST20  = []	
	gr_color_LIST20   = []

	RA_LIST19   = []      
	DEC_LIST19  = []		
	z_LIST19    = []		
	abs_g_mag_LIST19  = []	
	abs_r_mag_LIST19  = []	
	gr_color_LIST19   = []

	RA_LIST18   = []      
	DEC_LIST18  = []		
	z_LIST18    = []		
	abs_g_mag_LIST18  = []	
	abs_r_mag_LIST18  = []	
	gr_color_LIST18   = []

	for line in fp:
		current_index += 1
		if checkcondition(condition1, condition2) == False:
			print("Condition FALSE.")
			break
		else:
			pass		#Do nothing. Continue on with code. 


		splitline = line.split()  #Divide the non-hashed data into splitline.
		
		"""
		---------------------------------------
		# Pulls the fields from the data file. 
		---------------------------------------		
		"""
		RA_value   = float(splitline[0])        #Sets RA
		DEC_value  = float(splitline[1])		#Sets Dec
		z_value    = float(splitline[2])		#Sets  z
		abs_r_mag  = float(splitline[3])		#Sets  r-band mag
		abs_g_mag  = float(splitline[4])		#Sets  g-band mag

		"""
		---------------------------------------
		# Calculates values based on fields
		---------------------------------------
		"""
		g_r_color  = abs_g_mag - abs_r_mag
	
		"""
		---------------------------------------
		# Logic to count the number of blue and red galaxies. 
		---------------------------------------
		"""
		if g_r_color > 0.75:
			g_r_more7p5counter += 1
		else:
			g_r_less7p5counter += 1

		if z_value > 0. and abs_r_mag == -20.:
			print("Entering test loop: z > 0.15 and Mr=-20")


		"""
		---------------------------------------		
		#Logic that sets the index cutoff for the three subsample groups
		---------------------------------------
		"""
		if ((abs_r_mag <= -18.0) and (z_value <= z_lim18)):
			RA_LIST18.append(RA_value)
			DEC_LIST18.append(DEC_value)
			z_LIST18.append(z_value)
			abs_g_mag_LIST18.append(abs_g_mag)
			abs_r_mag_LIST18.append(abs_r_mag)
			gr_color_LIST18.append(g_r_color)
			#Counts the g-r for the limited sample
			if g_r_color >= 0.75:
				g_r_more7p5counter18 += 1
			else:
				g_r_less7p5counter18 += 1

		if ((abs_r_mag <= -19.0) and (z_value <= z_lim19)):
			RA_LIST19.append(RA_value)
			DEC_LIST19.append(DEC_value)
			z_LIST19.append(z_value)
			abs_g_mag_LIST19.append(abs_g_mag)
			abs_r_mag_LIST19.append(abs_r_mag)
			gr_color_LIST19.append(g_r_color)
			#Counts the g-r for the limited sample				
			if g_r_color >= 0.75:
				g_r_more7p5counter19 += 1
			else:
				g_r_less7p5counter19 += 1
		if ((abs_r_mag <= -20.0) and (z_value <= z_lim20)):
			#if (abs_r_mag == -20.00 or abs_r_mag == -20.01) and (z_value <= z_lim20):
			RA_LIST20.append(RA_value)
			DEC_LIST20.append(DEC_value)
			z_LIST20.append(z_value)
			abs_g_mag_LIST20.append(abs_g_mag)
			abs_r_mag_LIST20.append(abs_r_mag)
			gr_color_LIST20.append(g_r_color)
			#Counts the g-r for the limited sample
			if g_r_color >= 0.75:
				g_r_more7p5counter20 += 1
			else:
				g_r_less7p5counter20 += 1

		
		"""
		---------------------------------------
		APPEND VALUES TO TOTAL LISTS
		---------------------------------------
		"""
		RA_LIST.append(RA_value)
		DEC_LIST.append(DEC_value)
		z_LIST.append(z_value)
		abs_g_mag_LIST.append(abs_g_mag)
		abs_r_mag_LIST.append(abs_r_mag)
		gr_color_LIST.append(g_r_color)


"""
----------------------------------------------
# BEGIN PLOT FUNCTION Calls. { X,  Y}
----------------------------------------------
"""

#======================================================
# A.
#======================================================
title_label = "DEC vs. RA: FULL DR7"
x_label = "RA"
y_label = "DEC"
x_data  = RA_LIST
y_data  = DEC_LIST 
pointsize = 1
yflip = False

#pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)

#======================================================
# A. part II - M_r vs. Redshift
#======================================================
title_label = "r-band Mag vs. Redshift"
x_label = "Redshift, z"
y_label = "M_r"
x_data  = z_LIST
y_data  = abs_r_mag_LIST
pointsize = 1
yflip = True

#pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)
x_data[:] = []
y_data[:] = []
#======================================================
# B. Color Distribution of Galaxies
#======================================================
title_label = "g-r Color Distribution of galaxies"
x_label = "Redshift, z"
y_label = "g-r color"
x_data  = gr_color_LIST
y_data  = z_LIST
pointsize = 1
yflip = True

print("There were %s  BLUE galaxies!" % g_r_less7p5counter)
print("There were %s  RED  galaxies!" % g_r_more7p5counter)
#pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)
x_data[:] = []
y_data[:] = []
#======================================================
# C.  r-band luminosity function
#======================================================



#======================================================
# D.  Volume Limited sample { -20, -19, -18 }
#======================================================
#NOTE: THE FOLLOWING ASSUMES THAT THE DATA IS ORDERED ON M_r 
#-------------------------------------------
# Subsample -18
#-------------------------------------------
title_label = "r-band Mag vs. Redshift: -18"
x_label = "Redshift, z"
y_label = "M_r"
x_data  = z_LIST18
y_data  = abs_r_mag_LIST18
pointsize = 1
yflip = True
mag = -18
print_subsample_info(mag, x_data, g_r_more7p5counter18, g_r_less7p5counter18 )

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)

x_data[:] = []
y_data[:] = []

#-------------------------------------------
# Subsample -19
#-------------------------------------------
title_label = "r-band Mag vs. Redshift: -19"
x_label = "Redshift, z"
y_label = "M_r"
x_data  = z_LIST19
y_data  = abs_r_mag_LIST19
pointsize = 1
yflip = True
mag = -19
print_subsample_info(mag, x_data, g_r_more7p5counter19, g_r_less7p5counter19 )
pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)


x_data[:] = []
y_data[:] = []
#-------------------------------------------
# Subsample -20
#-------------------------------------------
title_label = "r-band Mag vs. Redshift: -20"
x_label = "Redshift, z"
y_label = "M_r"
x_data  = z_LIST20
#print abs_r_mag_LIST20
y_data  = abs_r_mag_LIST20
pointsize = 1
yflip = True
mag = -20
print_subsample_info(mag, x_data, g_r_more7p5counter20, g_r_less7p5counter20)
print("MAX z value of z_LIST20: %.4f " % max(z_LIST20) )
pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)


#print("M_r < -20:  \n ===============================\n")
#print("The Redshift bounds are: [ %.4f --> %.4f ] " % (min(x_data) , max(x_data)))
#print("The Volume is:             %.4f [(h^-1 Mpc)^3]." % get_volume(max(x_data)))
#print("The Red Blue Fraction is: %.3f" %  get_bluefr(g_r_more7p5counter20, g_r_less7p5counter20))

x_data[:] = []
y_data[:] = []
#======================================================
# E.  Luminosity Function of 3 Volume limited samples
#======================================================
"""
title_label = "Luminosity Function"
x_label = "x"
y_label = "y"
x_data  = gr_color_LIST
y_data  = z_LIST
pointsize = 1
yflip = True

print("There were %s  BLUE galaxies!" % g_r_less7p5counter)
print("There were %s  RED  galaxies!" % g_r_more7p5counter)
pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)

"""

"""
title_label = "TEST PLOT"
x_label = "x"
y_label = "y"
y_data  = abs_r_mag_LIST
x_data  = np.linspace(1,len(abs_r_mag_LIST), len(abs_r_mag_LIST))
pointsize = 1
yflip = True

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, pointsize, yflip, pcount)
"""



print("\nEND. No more things to do for the CPU!")







# Garbage comment for sublime convienence akjsdfjan


