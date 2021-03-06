import sys
import os
import numpy as np
import pandas as pd
import math
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

def maxmass_given_time(time_gyrs, init_maxmass, init_minmass):
	#Given a time since initial time, this gives what mass will be off-MS.
	#Time must be in solar units! Note: t/t_sun = (m / m_sun)**-3 
	if time_gyrs == 0:
		return init_maxmass

	t_sun    = 10.0  #Gyrs
	#max_mass_sol = (time_gyrs / t_sun) ** (-1./3.)  
	max_mass_sol = (time_gyrs / t_sun) ** (-1./2.5)  
	print("Max Mass: %.4f " % max_mass_sol)
	if max_mass_sol > init_maxmass:
		max_mass_sol = init_maxmass
	elif max_mass_sol < init_minmass:
		print("At current time, init_minmass has been reached. *Function: maxmass_given_time()")
	return max_mass_sol

def get_mag_given_luminosity(total_luminosity_at_time, M_bol_sun):
	#Calculates the absolute magnitude given a luminosity. 
	mag_bol_total = -2.5 * math.log10(total_luminosity_at_time) + M_bol_sun
	return mag_bol_total

def calc_mag_diff(mag_initial, mag_bol_total):
	mag_difference = mag_bol_total - mag_initial 
	return mag_difference

def get_g_minus_r_from_mass( m ):
	#Approximates (crudely) g-r given the mass of a star.  m in solar units. 
	g_minus_r = math.log((m+2)/m) - 0.65
	return g_minus_r 

def get_time_increment(time_increment, time):
	ti = time_increment*((time**(0.5)+0.1)/5.)
	#ti = time_increment
	return ti


#BEGIN MAIN PROGRAM
#-------------------------------------
print "Beginning Program!"
print "Initializing arrays and values ..."

#Define the steps the code will take: 
current_iteration = 0
iteration_steps   = 20000

iteration_masssteps = iteration_steps
iteration_timesteps = iteration_steps*0.1

#Time (units: Gyrs)
time_initial = 0 
time_final   = 11.

#Mass Range: (Units: Solar Mass Units)
mass_range_small = 0.08
mass_range_large = 50.0 

#Sets N_0 normalization. 
N_0 = 1.0

#Sets the Bol Mag of the sun. 
M_bol_sun = 4.745

#Sets the initial Mass Range for integration. Changes with time after initial step!
mass_initial   = mass_range_small
mass_final     = mass_range_large

mass_increment = calc_increment(mass_final, mass_initial, iteration_masssteps)
time_increment = calc_increment(time_final, time_initial, iteration_timesteps)

#Initialize Blank Arrays for plotting
mag_diff_time_array   = []
luminosity_time_array = []
time_array            = []
g_r_color_array       = []   #An array to accurately track the g-r color of IMF. 
g_r_color_est_array   = []   #An array to track the most massive star color

#Set t = 0 to start right after IMF. 
t = 0

#While Loop #2
while t <= time_final:
	#Sets the current luminosity to zero. Preparation for Loop over Mass range. 
	total_luminosity_at_time  = 0 
	g_r_weighted_by_lum       = 0 
	g_r_weighted_by_lum_total = 0
	m = mass_range_small  #This could be set to mass_initial but it shouldnt change. 
	
	#While Loop #1
	mass_final = maxmass_given_time(t, mass_range_large, mass_range_small)
	while m <= mass_final: 
		#Mass_initial is the small masses and mass_final are large masses.
		N_stars_massbin = N_0 * Salpeter_IMF(m, mass_increment)
		luminosity_given_mass = calc_Luminosity(m+0.5*mass_increment) 

		#Recursively add luminosities at mass bins to get a total luminosity at t = tj
		total_luminosity_at_time += (luminosity_given_mass * N_stars_massbin)

		g_r_given_mass = get_g_minus_r_from_mass( m )
		g_r_weighted_by_lum = g_r_given_mass * (luminosity_given_mass * N_stars_massbin)

		g_r_weighted_by_lum_total += g_r_weighted_by_lum

		#Update m for while loop 1
		m += mass_increment


	mag_bol_total = get_mag_given_luminosity(total_luminosity_at_time, M_bol_sun)
	if t == 0:
		mag_initial = mag_bol_total
		print("Initial Magnitude Set: %.3f " % mag_initial)

	#Set the Mass Difference to be plotted
	mag_difference = calc_mag_diff(mag_initial, mag_bol_total)
	mag_diff_time_array.append(mag_difference)

	#Set the g - r array to be plotted
	color_value = g_r_weighted_by_lum_total / total_luminosity_at_time
	g_r_color_est_array.append( g_r_given_mass )  #Appends the largest mass color at t. 
	g_r_color_array.append( color_value )

	luminosity_time_array.append(total_luminosity_at_time)
	time_array.append(t)

	#Update t for while loop 2
	time_increment2 = get_time_increment(time_increment, t) 
	t += time_increment2 



plt.figure(1)	
plt.plot(time_array, luminosity_time_array, linestyle='-', marker='')

plt.figure(2)
plt.plot(time_array, mag_diff_time_array, linestyle='-', marker='')

plt.figure(3)
plt.plot(time_array, g_r_color_array, linestyle='-', marker='', label='Average g-r')
plt.plot(time_array, g_r_color_est_array, linestyle='-', marker='', label='MaxMass g-r')

plt.figure(4)
plt.plot(g_r_color_array, mag_diff_time_array, linestyle='-', marker='', label='Average g-r')
plt.plot(g_r_color_est_array, mag_diff_time_array, linestyle='-', marker='', label='MaxMass g-r')


#Plot the relevant time steps. 
problem3_time_list   = [ 0.01, 0.1, 1., 2., 5., 10.]
p3_time_list_labels  = [ '10 Myrs', '100 Myrs', '1 Gyr', '2 Gyrs', '5 Gyrs', '10 Gyrs']
time_marker_list     = [ '>' , 'x', '^', 'o', '+', 'D']
time_markersize_list = [ 5, 7 , 9, 11, 13, 15]
#Loops to find and plot the relevant Ages. 
#for time_ in time_array:
j = 0
for k in range(0,len(problem3_time_list)):
	poop = problem3_time_list[k]
	print poop
	timepoop = time_array[j]  
	while poop >= timepoop:	
		j += 1
		timepoop = time_array[j]

	#if time_ == problem3_time_list[k]:
	plt.plot(g_r_color_array[j], mag_diff_time_array[j], linestyle='', \
			marker=time_marker_list[k], markersize=time_markersize_list[k],  label=str(p3_time_list_labels[k]) )
	plt.plot(g_r_color_est_array[j], mag_diff_time_array[j], linestyle='', \
			marker=time_marker_list[k], markersize=time_markersize_list[k] )



#PLOTTING
#-------------------------------------
figure_name="Jan31HW1P1LuminosityD.png"
plot_title="Luminosity vs. Time"  #Can code the number in with treemax
x_axis="Time [Gyrs]"
y_axis="Luminosity [solar units * N_0]"

figure_name2="Jan31HW1P1MagnitudesD.png"
plot_title2="Difference in Magnitude vs. Time"  #Can code the number in with treemax
x_axis2="Time [Gyrs]"
y_axis2="Absolute Magnitude Change: M(t) - M(t=0)"

figure_name3="Jan31HW1P2ColorsD.png"
plot_title3="Average Color vs. Time"  #Can code the number in with treemax
x_axis3="Time [Gys]"
y_axis3="Average Color: g - r"

figure_name4="Jan31HW1P3CompareColorMagD.png"
plot_title4="Mag Difference vs. Color"  #Can code the number in with treemax
y_axis4="Mag Difference M(t) - M(t=0)"
x_axis4="Average Color: g - r"


plt.figure(1)
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)
plt.xlim([0,10])
plt.yscale('log')
#plt.xscale('log')
plt.savefig(figure_name)

plt.figure(2)
plt.title(plot_title2)
plt.xlabel(x_axis2)
plt.ylabel(y_axis2)
plt.xlim([0,10])
#plt.xscale('log')
plt.savefig(figure_name2)

plt.figure(3)
plt.title(plot_title3)
plt.xlabel(x_axis3)
plt.ylabel(y_axis3)
plt.legend()
plt.xlim([0,10])
#plt.xscale('log')
plt.savefig(figure_name3)

plt.figure(4)
plt.title(plot_title4)
plt.xlabel(x_axis4)
plt.ylabel(y_axis4)
plt.legend(title='Ages', loc='best')
plt.savefig(figure_name4)


print("END.")

