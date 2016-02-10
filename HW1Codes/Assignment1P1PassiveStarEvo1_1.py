import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use( 'Agg' )
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

"""
Problem 1: Write a code to compute the passive evolution of a stellar population whose star
formation shuts off at time t=0. For the sake of simplicity, you can assume that all stars
are on the Main Sequence (i.e., when they die, they effectively disappear rather than
become giants)

Assume that at time t=0, the stellar population has a Salpeter mass function and contains
stars of all masses in the range 0.08 - 50 Msun ...
1. Calculate the total luminosity of the stellar population as a function of time in units of
its starting luminosity L0 . Then convert this to absolute magnitude: M(t)-M(0), the
difference between the absolute magnitude at time t and that at time t=0.
Plot M(t)-M(0) as a function of time (make the x-axis run from 0 to 10Gyr).

Steps: 
1. Define the Salpeter Mass Function.  IMF(m)*dm = IMF_init (m / M_sun)^-2.35 * (dm/M_sun)
	This states: The number of stars in each bin of mass is equal to the number of stars
	in the suns 'bin' of mass times the fractional mass to the minus 2.35, times the bin
	of mass normalized to the sun.  Note: IMF(m) = N(m) sometimes referred to by XI.

2. To get the total luminosity at initial time, integrate the IMF over all masses and compute.
	Note:  L / L_sun = (Mass / Mass_sun) ^ alpha,   where [1<alpha<6.]  ~3.5. 
	A variable alpha should be considered given that the range is outside of [2<mass<20].
	(The first generation of code will assume: alpha = 3.5.)
	

3. Recomute the bounds of integration given the main sequence lifetime, in bins of time. 
	Note: Ideally the bins of time would not be standardly spaced by would scale with
	the time since t=0 to account for the varying rate of decay. 
		* Assuming 0 --> 10 Gyrs, a spacing of 1 Myrs would lead to 10,000 bins. 
		
4. Integrate over the new bounds for the new luminosity until 10Gyr. 

5. Convert Luminosity to absolute Magnitude and plot Results. 
	M_bol = -2.5 * log_b10[( L / L_sun )]. 

6. Plot the Magnitude difference as a function of time.  
"""


#Function:
"""
#FUNCTIONS
#-------------------------------------
"""

def calc_Luminosity(mass):
	#Recieves a mass in solar units and returns the Luminosity in solar units.
	alpha = 3.5
	lum_total_sun = (mass)**alpha
	return lum_total_sun


def Salpeter_IMF(mass, dm):
	#This function computes the Salpeter IMF value, given a mass and (dm). [sol_units]
	alpha = -2.35
	IMF_value = ((mass)**(alpha))*dm
	return IMF_value

def calc_increment(big_mass, small_mass, increment):
	#Calculates the difference in solar units in a linear fashion.
	mass_difference = (big_mass - small_mass)/increment
	return mass_difference

#BEGIN MAIN PROGRAM
#-------------------------------------
print "Beginning Program!"
print "Initializing arrays and values ..."

#Define the steps the code will take: 
current_iteration = 0
iteration_steps   = 1000

iteration_masssteps = iteration_steps
iteration_timesteps = iteration_steps

#Time (units: Gyrs)
time_initial = 0 
time_final   = 10.0

#Mass Range: (Units: Solar Mass Units)
mass_range_small = 0.08
mass_range_large = 20.0 

#Sets N_0 normalization. 
N_0 = 1.0

#Sets the initial Mass Range for integration. Changes with time after initial step!
mass_initial   = mass_range_small
mass_final     = mass_range_large

mass_increment = calc_increment(mass_final, mass_initial, iteration_masssteps)
time_increment = calc_increment(time_final, time_initial, iteration_timesteps)


luminosity_time_array = []
time_array = []

#While Loop #2
t = 0
while t <= time_final:

	#Sets the current luminosity to zero. Preparation for Loop over Mass range. 
	total_luminosity_at_time = 0  
	m = mass_range_small   #This could be set to mass_initial but it shouldnt change. 
	#While Loop #1
	while m <= mass_final: 
		#Mass_initial is the small masses and mass_final are large masses.
		N_stars_massbin = N_0 * Salpeter_IMF(m, mass_increment)
		luminosity_given_mass = calc_Luminosity(m+0.5*mass_increment) 

		#Recursively add luminosities at mass bins to get a total luminosity at t = tj
		total_luminosity_at_time += (luminosity_given_mass * N_stars_massbin)
		#Update m for while loop 1
		m += mass_increment

	luminosity_time_array.append(total_luminosity_at_time)
	time_array.append(t)
	#Update t for while loop 2
	t += time_increment 



plt.figure(1)	
plt.plot(time_array, luminosity_time_array, linestyle='-', marker='o')


#plt.figure(2)
#plt.plot(r_arr, press_arr, linestyle='-', marker='')

#PLOTTING
#-------------------------------------
figure_name="HW1P1testfigure.png"
plot_title="null"  #Can code the number in with treemax
x_axis="null"
y_axis="loops"

"""
figure_name2="figure2hw2.png"
plot_title2="Pressure vs. r"  #Can code the number in with treemax
x_axis2="r"
y_axis2="Pressure"

figure_name3="figure1hw3.png"
plot_title3="Theta' vs. Xi"  #Can code the number in with treemax
x_axis3="Xi"
y_axis3="Theta''"
"""

plt.figure(1)
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)
plt.savefig(figure_name)

"""
plt.figure(2)
plt.title(plot_title2)
plt.xlabel(x_axis2)
plt.ylabel(y_axis2)
plt.savefig(figure_name2)
"""
print("END.")

