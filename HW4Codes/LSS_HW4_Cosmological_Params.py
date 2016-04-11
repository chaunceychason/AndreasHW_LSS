import sys
import os
import numpy as np
#import pickle
#import pandas as pd
#from scipy import stats  
#import astropy.units as u 
from astropy import cosmology
#import astropy 
from astropy.cosmology import w0waCDM 

import matplotlib
matplotlib.use( 'Agg' )
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

"""
Program Notes:

HW AST8050: Computing Cosmological Quantities. 

FOR FLAT UNIVERSES ONLY! 

* Input {hubble constant, matter density, dark energy density, dark energy eq. of state}
        {h, Om, Oa, w}
* Output {covmoving distance, Luminosity distance, 
 ang. diamter distance, comoving volume, lookback time, - vs. redshift  }

Notes: 
- Redshift Range [0,3]

- distnace units h^-1 Gpc.  volume (d^3)

- Lookback time in Gyrs. 
"""

#FUNCTIONS
#--------------------------
def getvol(radius):
	nothing = 0
	vol = 4./3. * np.pi * radius**3
	return vol


#BEGIN MAIN PROGRAM
#-------------------------------------
print "Beginning HW4 LSS Program!"
print "This code will determine cosmological parameters"
print
print "Setting Date String"
datestring = "Apr11"

num_z_points = 100
redshift_max = 3
redshift_min = 0
redshifts = np.linspace(redshift_min, redshift_max, num=num_z_points)


# h, Om, Oa, w 
param1 = np.array([1,  1.0, 0.0,   0  ])
param2 = np.array([1, 0.25, 0.75, -1  ])
param3 = np.array([1, 0.25, 0.75, -0.8])
param4 = np.array([1, 0.25, 0.75, -1.2])
hval = 1.0

paramALL = np.array([param1, param2, param3, param4])

cosmo1 = w0waCDM(H0=100.*hval*param1[0], Om0=param1[1], Ode0=param1[2] )
cosmo2 = w0waCDM(H0=100.*hval*param2[0], Om0=param2[1], Ode0=param2[2], w0=param2[3])
cosmo3 = w0waCDM(H0=100.*hval*param3[0], Om0=param3[1], Ode0=param3[2], w0=param3[3])
cosmo4 = w0waCDM(H0=100.*hval*param4[0], Om0=param4[1], Ode0=param4[2], w0=param4[3])

"""
test_z = 0.5
dc = cosmo1.comoving_distance(test_z)
print("Test value for dc = ", dc)
"""
#Units: Mpc
dcmpc1 = cosmo1.comoving_distance(redshifts)
dcmpc2 = cosmo2.comoving_distance(redshifts)
dcmpc3 = cosmo3.comoving_distance(redshifts)
dcmpc4 = cosmo4.comoving_distance(redshifts)

#Units: Mpc
ldmpc1 = cosmo1.luminosity_distance(redshifts)
ldmpc2 = cosmo2.luminosity_distance(redshifts)
ldmpc3 = cosmo3.luminosity_distance(redshifts)
ldmpc4 = cosmo4.luminosity_distance(redshifts)

#Units: Gyr
age1 = cosmo1.age(redshifts) 
age2 = cosmo2.age(redshifts) 
age3 = cosmo3.age(redshifts) 
age4 = cosmo4.age(redshifts) 

#Units: Gyrs
lbt1 = cosmo1.lookback_time(redshifts) 
lbt2 = cosmo2.lookback_time(redshifts) 
lbt3 = cosmo3.lookback_time(redshifts) 
lbt4 = cosmo4.lookback_time(redshifts) 

#Units: Mpc
admpc1 = cosmo1.angular_diameter_distance(redshifts)
admpc2 = cosmo2.angular_diameter_distance(redshifts)
admpc3 = cosmo3.angular_diameter_distance(redshifts)
admpc4 = cosmo4.angular_diameter_distance(redshifts)


#CONVERT TO APPROPIATE UNITS -> Mpc to Gpc.

dc1 = dcmpc1 / 1000.
dc2 = dcmpc2 / 1000.
dc3 = dcmpc3 / 1000.
dc4 = dcmpc4 / 1000.

ld1 = ldmpc1 / 1000.
ld2 = ldmpc2 / 1000.
ld3 = ldmpc3 / 1000.
ld4 = ldmpc4 / 1000.

ad1 = admpc1 / 1000.
ad2 = admpc2 / 1000.
ad3 = admpc3 / 1000.
ad4 = admpc4 / 1000.


#COMOVING PLOT
plot_title="Comoving Distance vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Comoving Distance Gpc/h"
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(redshifts, dc1, label="cosmo1")
plt.plot(redshifts, dc2, label="cosmo2")
plt.plot(redshifts, dc3, label="cosmo3")
plt.plot(redshifts, dc4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4dc' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()

#LUMINOSITY DISTANCE PLOT
plot_title="Luminosity Distance vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Luminosity Distance Gpc/h"
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(redshifts, ld1, label="cosmo1")
plt.plot(redshifts, ld2, label="cosmo2")
plt.plot(redshifts, ld3, label="cosmo3")
plt.plot(redshifts, ld4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4ld' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()

#ANG DIAMETER DISTANCE PLOT
plot_title="Ang. Diameter Distance vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Ang. Diameter Distance "
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(redshifts, ad1, label="cosmo1")
plt.plot(redshifts, ad2, label="cosmo2")
plt.plot(redshifts, ad3, label="cosmo3")
plt.plot(redshifts, ad4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4ad' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()


#COMOVING VOLUME PLOT
plot_title="Comoving Volume vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Volume Gpc^3/h^3"
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

vol1 = getvol(dc1)
vol2 = getvol(dc2)
vol3 = getvol(dc3)
vol4 = getvol(dc4)
plt.plot(redshifts, vol1, label="cosmo1")
plt.plot(redshifts, vol2, label="cosmo2")
plt.plot(redshifts, vol3, label="cosmo3")
plt.plot(redshifts, vol4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4vol' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()


#AGE DISTANCE PLOT
plot_title="Age vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Age Gyrs"
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(redshifts, age1, label="cosmo1")
plt.plot(redshifts, age2, label="cosmo2")
plt.plot(redshifts, age3, label="cosmo3")
plt.plot(redshifts, age4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4age' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()


#Lookback Time PLOT
plot_title="Lookback Time vs. Z"   #Can code the number in with treemax
x_axis="z (redshift)"
y_axis="Lookback Time Gyrs/h"
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

plt.plot(redshifts, lbt1, label="cosmo1")
plt.plot(redshifts, lbt2, label="cosmo2")
plt.plot(redshifts, lbt3, label="cosmo3")
plt.plot(redshifts, lbt4, label="cosmo4")
figure_name=os.path.expanduser('~/LSSHW4lbt' +'.png')
plt.legend(loc='best')
plt.savefig(figure_name)
print("Saving plot: %s" % figure_name)
plt.clf()





print("Program Completed. ")
