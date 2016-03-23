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
#	    Galaxy Samples:  LENGTHS = { 5495, 28162, 28383}
#	    { SDSS_Mr21_rspace.dat, SDSS_Mr20_rspace.dat, SDSS_Mr20_zspace.dat }
#	    Random Points: 
#	    { SDSS_random.dat }  LENGTH = 42654
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

r_min = 0.1 
r_max = 20
logr_min = -1
logr_max = 1.301


datafilename1 = 'output_Xi_DM.txt'
datafilename2 = 'output_Xi_20r.txt'
datafilename3 = 'output_Xi_21r.txt'


#with open('time_stats_run1temp.txt' , "r") as f2:

r_list = np.logspace(logr_min, logr_max, 15)
logr_list = np.log10(r_list)
logr_datacountlist = np.zeros(15) 

#CHECK IF FILE EXISTS.
if os.path.exists(datafilename1):
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")    

if os.path.exists(datafilename2):
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")  
if os.path.exists(datafilename3):
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")  

#Read in datafile:
print("Reading File: Data should be structured...")
Xi1data = np.loadtxt(datafilename1)
Xi2data = np.loadtxt(datafilename2)
Xi3data = np.loadtxt(datafilename3)

N_sample = len(Xi1data)
print("Lenth of XiData is %d" % N_sample)


print("Made the Arrays!")
print("======================")
print("For values of logr of:")
print logr_list 
print("For values of r of:")
print r_list 

print("THE LOGDATACOUNT is:")
print logr_datacountlist 

"""
theoutfile = open('random_logr_datacounts.txt', 'w')
for item in logr_datacountlist:
	print>>theoutfile, item
theoutfile.close()
"""


#COMPUTE BIAS!

bias1 = np.sqrt(Xi2data/Xi1data)  #20r
bias2 = np.sqrt(Xi3data/Xi1data)  #21r


#THIS PROGRAM PLOTS DATA IN FORM 'x, y'
"""
=====================================
    Correlation Function Plot with DM
=====================================
"""
plot_title="Correlation Function -20 & -21"
x_axis="log Distance [h-1 Mpc]"
y_axis="Xi(r)"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)


plt.plot(logr_list, Xi1data, color='b' ,marker='o', label='Mr > -20r' )
#plt.plot(logr_list, Xi3data, color='g' ,marker='o', label='Mr > -21r' )

#plt.xscale('log')
#plt.yscale('log')
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HW3CorrelationFuncParallelDM.png"
plt.savefig(tmp_filename, rasterized=True)
plt.clf()

"""
=========================================
    BIAS Function Plot 
=========================================
"""

plot_title="Bias for -20 real vs z"
x_axis="Log Distance [h-1 Mpc]"
y_axis="Bias (r)"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)


plt.plot(logr_list, bias1, color='b' ,marker='o', label='Mr > -20r' )
plt.plot(logr_list, bias2, color='g' ,marker='o', label='Mr > -21r' )

#plt.xscale('log')
#plt.yscale('log')
#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='best')

#Saves Plot
tmp_filename = "HW3CorrelationFuncParallelBias.png"
plt.savefig(tmp_filename, rasterized=True)
plt.clf()

print "The program has finished running. All files closed. \nThe results should be in your directory"



#End