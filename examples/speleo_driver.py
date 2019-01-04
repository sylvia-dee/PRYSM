#!/usr/bin/env python
# This is a driver script for the Proxy System Model
# of the d18O of  Speleothem Calcite
# Created 10/28/2014 by Julien Emile-Geay (julieneg@usc.edu)
# Modified 01/27/2015 by Sylvia Dee <sdee@usc.edu>
# Modified 11/17/2015 by Sylvia Dee <sylvia_dee@brown.edu>
#==============================================================
"""
Driver for Speleothem Calcite  Model
~ Load Environmental Variables to compute d18O calcite/dripwater.

Required Isotope-Enabled GCM INPUTS for PSM:

d18O (Soil Water)
d18O (Precipitation)
Precipitation (mm/day)
Temperature (K)

"""
#==============================================================
# E0. Initialization
#==============================================================
import numpy as np
import matplotlib.pyplot as plt
import psm.aux_functions.butter_lowpass_filter as bwf
#import scipy.interpolate as si
from psm.speleo.sensor import * #We should probably be explicit
from psm.aux_functions.analytical_error import analytical_error
from psm.aux_functions.analytical_err_simple import analytical_err_simple
from psm.agemodels.tiepoint import chron
#R things. It would be nice to get rid of this dependency
import rpy2.robjects as robjects
from rpy2.robjects import FloatVector
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

#==============================================================
# E1. Load input data
#==============================================================
print('Loading data...')
# Note: for the tester files, the data corresponds to Borneo
# output is from SPEEDY-IER
# USER should please load your own data here!
datadir='test_data_speleo/'
cave_data = {}
cave_data['dO18_prcp'] = np.load(datadir+'B_PISO.npy')
cave_data['dO18_soil'] = np.load(datadir+'B_SISO.npy')
cave_data['temp']      = np.load(datadir+'B_TEMP.npy')
cave_data['prcp']      = np.load(datadir+'B_PRCP.npy')

# define time axis, and the time stepping parameter dt
t=np.arange(1000,2005,1.0)
# define time step in years. for annual data, dt = 1.0. for monthly data, dt=1./12.
dt=1.0
#dt=1./12.
#==========================================================================
# E2. APPLY SENSOR MODEL
#==========================================================================
print('Running sensor model...')
model = 'Adv-Disp'
Pe   = 1.0  # Peclet number
tau0 = 1.0  # mean aquifer transit time
d18Op = np.array(cave_data['dO18_prcp'])
d18Os = np.array(cave_data['dO18_soil'])
T     = np.array(cave_data['temp'])
P     = np.array(cave_data['prcp'])

#  Apply two different catchment models with Tau = 1.0 years, 5.0 years
d18O_wm1, d18O_wmK1, hwm1 = speleo_sensor(dt,d18Os,T,'Well-Mixed',1.0)
d18O_ad1, d18O_adK1, had1 = speleo_sensor(dt,d18Os,T,'Adv-Disp',1.0,Pe)
d18O_wm2, d18O_wmK2, hwm2 = speleo_sensor(dt,d18Os,T,'Well-Mixed',5.0)
d18O_ad2, d18O_adK2, had2 = speleo_sensor(dt,d18Os,T,'Adv-Disp',5.0,Pe)

# Save a Zipped output file (optional)
# >np.savez('cave_dataCave_speleothem_example',cave_data, d18O_c1, d18O_K1,d18O_c2, d18O_K2)

# lowpass filter everything using a zero-phase butterworth filter
fc = 0.1 # cutoff frequency
Tl = bwf.filter(T,fc)
Pl = bwf.filter(P,fc)
d18Opl = bwf.filter(d18Op,fc)
d18Osl = bwf.filter(d18Os,fc)
K1     = bwf.filter(d18O_wm1,fc)
K2     = bwf.filter(d18O_ad1,fc)
K3     = bwf.filter(d18O_wm2,fc)
K4     = bwf.filter(d18O_ad2,fc)

#==========================================================================
# E3. APPLY OBSERVATION MODEL
#==========================================================================
print('Running observation model...')

r = robjects.r

# apply analytical errors
#d18O_a = analytical_error(d18O_ad2,0.1,10)

# apply chronological errors
# data from Rasmussen et al, doi:10.1029/2006GL025714, 2006
d=60.# total length of core in cm
nyears = t.shape[0]
ages = FloatVector([-55,55,236,432,835])# age estimate
sd  = FloatVector([1,22,50,13,25])# SD of ages
positions = FloatVector([0,5,15,35,55])# position in core in cm
calCurves = StrVector(['normal','normal','normal','normal','normal'])
predictPositions = r.seq(0,d,by=d/nyears) # note d= total depth of core (60cm)

# Specify extractDate (top-most age of the core, in years BP)
topDate= -55

# Call BCHRON observation model to return CHRONS (in years BP)

chronBP, depth_horizons = chron(ages,sd,positions,calCurves,predictPositions,topDate)

# recast in years CE
chronCE = np.fliplr(1950 - chronBP)

#==========================================================================

# 3.2: Analytical Uncertainty Model:
print('Adding uncertainty...')
#Enter final simulated speleothem record (choose from above options)
X=d18O_wm1
#X=d18O_wm2

#3.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
speleo_upper, speleo_lower = analytical_err_simple(X,sigma)

#3.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
speleo_Xn=analytical_error(X,sigma)

#====================================================================
# Save whatever needs to be saved
print('Saving data...')
outdir='./results/'
np.save(outdir+"speleo_Xn.npy",speleo_Xn)
#Whatever else...
#====================================================================

###Move this to a separate plotting script ###
#f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
#f.set_figheight(12); f.set_figwidth(8)
#plt.rcParams['text.usetex']=True
#plt.rcParams['text.latex.unicode']=True
#plt.rcParams['font.family']='serif'
#
##  PHYSICAL VARIABLES
#ax1.spines["top"].set_visible(False)
#ax1.spines["right"].set_visible(False)
## Ensure that the axis ticks only show up on the bottom and left of the plot.
## Ticks on the right and top of the plot are generally unnecessary chartjunk.
#ax1.get_xaxis().tick_bottom()
#ax1.get_yaxis().tick_left()
#plt.yticks(fontsize=10); plt.xticks(fontsize=10)
## Remove the tick marks; they are unnecessary with the tick lines we just plotted.
#ax1.tick_params(axis="both", which="both", bottom="on", top="off",
#                labelbottom="on", left="on", right="off", labelleft="on",direction="out")
#ax1.grid(axis='y')
#
#ax1.plot(t,Tl-273.15,lw=1, color='Chocolate', label='Temp')
## Make the y-axis label and tick labels match the line color.
#ax1.set_ylabel('exp', color='Chocolate')
#for tl in ax1.get_yticklabels():
#    tl.set_color('Chocolate')
#ax1.set_xlim(1000,2000); ax1.set_ylim(13,21)
## label axes
#ax1.set_ylabel(r'Temp (${}^\circ$C)',fontsize=11)
#ax1.set_title('Climate Parameters',fontsize=14)
#
## add precip !!
#ax1b = ax1.twinx()
#ax1b.plot(t,Pl,lw=1, color='Teal')
#for tl in ax1b.get_yticklabels():
#    tl.set_color('Teal')
#ax1b.set_xlim(1000,2000); ax1b.set_ylim(0,1.5)
## label axes
#ax1b.set_ylabel('Precip (mm/day)',fontsize=11,color='Teal')
#
## GEOCHEMISTRY
#ax2.spines["top"].set_visible(False)
#ax2.spines["right"].set_visible(False)
#ax2.get_xaxis().tick_bottom()
#ax2.get_yaxis().tick_left()
#
## Remove the tick marks; they are unnecessary with the tick lines we just plotted.
#ax2.tick_params(axis="both", which="both", bottom="on", top="off",
#                labelbottom="on", left="on", right="off", labelleft="on",direction="out")
#ax2.grid(axis='y')
#
#ax2.plot(t,d18Opl,lw=1, color='Navy', label='Precip')
#ax2.plot(t,d18Osl,lw=1, color='DodgerBlue', label='Soil')
#
#lg = ax2.legend(fontsize=10,loc='center'); lg.draw_frame(False)
#ax2.set_xlim(1000,2000)
## label axes
#ax2.set_ylabel(r'$\delta^{18}O_\mathrm{water}$ (\u2030)',fontsize=12)
#ax2.set_title('Oxygen isotope source signals',fontsize=14)
#
## SENSOR MODELs
#ax3.spines["top"].set_visible(False)
#ax3.spines["right"].set_visible(False)
#ax3.get_xaxis().tick_bottom()
#ax3.get_yaxis().tick_left()
#
## Remove the tick marks; they are unnecessary with the tick lines we just plotted.
#ax3.tick_params(axis="both", which="both", bottom="on", top="off",
#                labelbottom="on", left="on", right="off", labelleft="on",direction="out")
#ax3.grid(axis='y')
#
#ax3.plot(t,K1,lw=1, color='DarkGreen', label=r'Well-mixed, $\tau_0 = 1$y')
#ax3.plot(t,K2,lw=2, linestyle = '-.', color='DarkGreen', label=r'Well-mixed, $\tau_0 = 5$y')
#ax3.plot(t,K3,lw=1, color='LimeGreen', label=r'Advection-Dispersion, $P_e = 1.0$, $\tau_0 = 1$y')
#ax3.plot(t,K4,lw=2, linestyle = '-.', color='LimeGreen', label=r'Well-mixed, $\tau_0 = 5$y')
#lg = ax3.legend(fontsize=10); lg.draw_frame(False)
#
#
#ax3.set_xlim(1000,2000); ax3.set_ylim(-9.5,-3.5)
## label axes
#ax3.set_xlabel('Time (year CE)',fontsize=12)
#ax3.set_ylabel(r'$\delta^{18}O_\mathrm{calcite}$ (\u2030)',fontsize=12)
#ax3.set_title('Karst Effects',fontsize=14)
#
#f.savefig('sensor_model_demo.pdf', dpi=400, facecolor='w', edgecolor='w',transparent=True)
#
#
#  relationship between d18Oc and precip

#f2, ax = plt.subplots(2, 2)
#f2.set_figheight(6); f2.set_figwidth(8)
# d18O rainfall vs rainfall
#ax[0,0].scatter(Pl,d18Opl,color='IndianRed',alpha=0.5)
#ax[0,0].set_xlabel('precipitation rate (mm/day)')
#ax[0,0].set_ylabel(r'$\delta^{18}\mathrm{O}_\mathrm{P}$ (\u2030)')
#plt.show()
#plt.title('Signal transformation on decadal time scales')
#f2.savefig('signal_transformation.pdf', dpi=200, facecolor='w', edgecolor='w',transparent=True)
