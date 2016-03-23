import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

"""
#  .-----..-.  .-. .---.              .-.-.  .---. .-..-. .-..-----.    
#  `-' '-'| {  } |/ {-. \     ___     | } }}/ {-. \{ ||  \{ |`-' '-'    
#    } {  {  /\  }\ '-} /    {___}    | |-' \ '-} /| }| }\  {  } {      
#    `-'  `-'  `-' `---'              `-'    `---' `-'`-' `-'  `-'      
#  .----. .---. .---. .---. .----..-.     .--. .-----..-. .---. .-. .-. 
#  | }`-'/ {-. \} }}_}} }}_}} |__}} |    / {} \`-' '-'{ |/ {-. \|  \{ | 
#  | },-.\ '-} /| } \ | } \ } '__}} '--./  /\  \ } {  | }\ '-} /| }\  { 
#  `----' `---' `-'-' `-'-' `----'`----'`-'  `-' `-'  `-' `---' `-' `-' 
#  .----..-. .-..-. .-..----..-----..-. .---. .-. .-.                   
#  } |__}| } { ||  \{ || }`-'`-' '-'{ |/ {-. \|  \{ |                   
#  } '_} \ `-' /| }\  {| },-.  } {  | }\ '-} /| }\  {                   
#  `--'   `---' `-' `-'`----'  `-'  `-' `---' `-' `-'                   
#                                                                       
# This program will make the correlation function given files with [RA, Dec, Z]. 
# It will use the estimator:
#      {   Xi(r) = (N_r/N_d)**2 (DD(r)/RR(r) - 1)   }
#
# The Sample data is separated into 15 logarithmic bins of separation in range 
# 0.1 - 20 h**-1 Mpc ( the first bin at log r = -1 and the last bin is at 
# log r = 1.301 )
#
# The Random Points are redshifts in range 0.02 <= z <= 0.06. Using the SDSS sky
# coverage. 
#
#
# FILE NAMES: 
#    Part 1:
#	    Galaxy Samples: 
#	    { SDSS_Mr21_rspace.dat, SDSS_Mr20_rspace.dat, SDSS_Mr20_zspace.dat }
#	    Random Points: 
#	    { SDSS_random.dat }
#
#		A.  log Xi vs. log r  for two real-space galaxies. {-20, -21}
#		B.  log Xi vs. log r  for real-space vs. redshift-space sample. {-20}
#		
#	Part 2: Use the DM files to find Xi(r).
#		DM Samples: Points in Cartesian Coordinates {Xi, Yi, Zi,} Vol=141.3h-1Mpc
#		{DM.dat}
#		Random Points:
#		{DM_random.dat}
#		
#		B. Compute the bias function: b(r) = SQRT[ Xi_gal / Xi_DM  ]
"""

def calc_Xi_r(  ):
	Xi = 1

	return Xi 

def radius_given_z(z):
	radius = 2950.*z
	return radius

def XYZ_given_RADECZ(ra, dec, z):
	#Returns Radius in {X, Y, Z} format in units of [h-1 Mpc]
	r = radius_given_z(z)
	PI_const = 3.14159
	dec_to_rad = ( 2.*np.pi/ 360.  )
	DEC_rad = dec * ( dec_to_rad )
	RA_rad  = ra  * ( dec_to_rad )
	x = r * math.cos( DEC_rad ) * math.cos( RA_rad)
	y = r * math.cos( DEC_rad ) * math.sin( RA_rad)
	z = r * math.sin( DEC_rad )
	return x, y, z

def distance_given_2points(x1, y1, z1, x2, y2, z2):
	#Distance in any units given. Such as Mpc.
	distance = ( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )**(1./2.)
	return distance

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
	#radius = 2950.*max_z 
	radius = radius_given_z(max_z)
	DR7_SAMPLE_COVERAGE_STERADIANS = 2.295   #Should be set to 2.295 sterad. (7675.2 deg^2)
	fract_vol = DR7_SAMPLE_COVERAGE_STERADIANS / (4.*np.pi)    #Fractional sky coverage DR7 to 4pi. 
	vol = fract_vol*(4./3.) * np.pi * (radius**3.) 
	return vol  #Units: [h-1 Mpc^3]

#Set to -1 and 1.301 as per instructions. 
logr_min = -1
logr_max = 1.301
errorcounts = 0

datafilename = 'SDSS_random.dat'
#with open('time_stats_run1temp.txt' , "r") as f2:

#logr_max = 1.301
print("Setting logr_max to a value NOT specificed in HW assignment. FIX!")

r_list = np.logspace(logr_min, logr_max, 15)
logr_list = np.log10(r_list)
logr_datacountlist = np.zeros(15) 
#CHECK IF FILE EXISTS.
if os.path.exists(datafilename):
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")    


#Read in datafile:
print("Reading File: Data should be structured (RA, DEC, Z)...")
RADECZdata = np.loadtxt(datafilename)

N_sample = len(RADECZdata)
#N_random = len(RandomData)

#Initialize X and Y arrays. 
xdata  = []
ydata1 = []

#//Initialize Distances Array.
#//distances = [0 for i in range(0, len(RADECZdata))]

print_5bool = 0
print_10bool = 0
print_20bool = 0
print_30bool = 0
print_40bool = 0
print_50bool = 0
print_75bool = 0
#print_90bool = 0
"""
=====================================
 CALCULATE THE CORRELATION FUNCTION
=====================================
"""
print("The length of RADECZdata is %d" % len(RADECZdata))
for i in range(0, len(RADECZdata)):
	x1,y1,z1 = XYZ_given_RADECZ(RADECZdata[i][0], RADECZdata[i][1], \
		RADECZdata[i][2])

	if i > 0.05*len(RADECZdata) and print_5bool == 0:
		print("Finished with 5 percent of the outer i loop. ")
		print_5bool = 1

	if i > 0.10*len(RADECZdata) and print_10bool == 0:
		print("Finished with 10 percent of the outer i loop. ")
		print_10bool = 1

	if i > 0.20*len(RADECZdata) and print_20bool == 0:
		print("Finished with 20 percent of the outer i loop. ")
		print_20bool = 1

	if i > 0.50*len(RADECZdata) and print_50bool == 0:
		print("Finished with 50 percent of the outer i loop. ")
		print_50bool = 1

	if i > 0.75*len(RADECZdata) and print_75bool == 0:
		print("Finished with 75 percent of the outer i loop. ")
		print_75bool = 1

	#print x1, y1, z1 
	for j in range(0, len(RADECZdata)):
		if i != j:
			x2,y2,z2 = XYZ_given_RADECZ(RADECZdata[j][0], RADECZdata[j][1], \
				RADECZdata[j][2]) 
			D = distance_given_2points(x1, y1, z1, x2, y2, z2)
			#print D
			logD = math.log10(D)
			#print logD
			#CHECK THE DISTANCE VALUE
			"""
			if -1.1 <= logD <= 1.303:
				pass
			else:
				print("It appears the distance of logD: %.3e is out of range" % logD)
			"""
			try:
				r_idx = np.min(np.where(logr_list > logD))
				logr_datacountlist[r_idx] += 1
			except:
				errorcounts += 1


print("FOUND THE COUNTS! Error Counts %d" % errorcounts)
print("======================")
print("For values of logr of:")
print logr_list 
print("For values of r of:")
print r_list 

print("THE LOGDATACOUNT is:")
print logr_datacountlist 

theoutfile = open('random_logr_datacounts.txt', 'w')
for item in logr_datacountlist:
	print>>theoutfile, item
theoutfile.close()

#THIS PROGRAM PLOTS DATA IN FORM 'x, y'
"""
=====================================
    Correlation Function Plot 
=====================================
"""


plot_title="xCorrelation Functionx "
x_axis="Log(Distance) [h-1 Mpc]"
y_axis="Xi(r)"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(logr_list, logr_datacountlist, color='b' ,marker='o', label='Random' )
plt.xscale('log')
plt.yscale('log')
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HW3CorrelationFuncRand.png"
plt.savefig(tmp_filename, rasterized=True)
plt.clf()





print "The program has finished running. All files closed. \nThe results should be in your directory"



#End