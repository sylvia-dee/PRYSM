#======================================================================
# Sylvia Dee
# Modified 03/08/16 <sylvia_dee@brown.edu>
# Modified 04/02/18 <sylvia@ig.utexas.edu>
# Script to run lake sediment proxy system model
#======================================================================

# PRYSM v2.0 Lake Sediments, DRIVER SCRIPT
# Notes: please see all function docstrings for proper use.
# Dependencies: rpy2, f2py, numpy

"""
    Lake PSM ~ Load Environmental Variables to compute lake energy/water balance.
    This is a driver script. It is set up to walk the user through the necessary
    commands to use the lake PSM. The user can load there own data in the fields below,
    or use the provided test data.

    This driver script calls the submodels:
        lake_sensor
        lake_archive
        lake_obs

    The final output is a pseudoproxy time series at a given site, with 1000
    plausible chronologies, age-depth information, and related sensor uncertaities.

"""
#======================================================================
# L0. Initialization
#======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from math import pi, sqrt, exp


#==============================================================
# L1. ENVIRONMENT SUB-MODEL:  Load environment Input data
#==============================================================

# cd /disk/staff/sylvia/LAKEPSM/

# # Example data directory
# datadir='test_data_lake/'
# print 'Loading data from ', datadir,' ...'

# env_input='example_input.txt'

#==============================================================

# F2PY to run environment model: run the following in the command window
# NOTE:
# (1) IT RUNS IN THE COMMAND LINE NOT IN THE PYTHON SHELL.
# (2) Make sure to run and install in the correct environment.

# f2py -c -m lakepsm lake_environment.f90

# MUST BE IN SAME DIRECTORY AS COMPILED FORTRAN WRAPPER

import lakepsm

lake_data_file= env_input

# Run Environment Model
lakepsm.lakemodel()

# The above should run the fortran model through a wrapper (F2Py)
# The output files (listed below) should now appear in your working directory.

# output files:
#profile_output.dat
#surface_output.dat

#==============================================================
# L2. SENSOR MODEL(S)
#==============================================================
import sensor_gdgt as gdgt
import sensor_carbonate as carb
import sensor_leafwax as leafwax

# Pull data (lake surface temperature) from output from env
# model and use to run sensor model
surf_ouput = 'ERA-HIST-Tlake_surf.dat'
surf_tempr = []
with open(surf_ouput, 'r') as data:
	tempr_yr = [] 
	for line in data:
		line_vals = line.split()
		tempr_yr.append(line_vals[1])
	surf_tempr.append(tempr_yr[341:-11])
surf_tempr = np.array(surf_tempr[0], dtype = float) +275.15

climate_input = 'ERA_INTERIM_climatology_Tanganyika.txt'
air_tempr = []
with open(climate_input, 'r') as data:
    airtemp_yr = []
    for line in data:
        line_vals = line.split()
        airtemp_yr.append(line_vals[2])
    air_tempr.append(airtemp_yr[:-12])
air_tempr = np.array(air_tempr[0], dtype = float)


LST = surf_tempr  # THIS SHOULD BE THE FULL TIME SERIES of lake surface temp
MAAT = air_tempr  # Full timeseries of air temperatures
beta = 1./50.
d18Ow = -2 # To be set by the user if isoflag is off

print('Running sensor model...')

# 2.1 GDGT Sensors
gdgt_proxy = gdgt.gdgt_sensor(LST,MAAT,beta=beta,model = 'MBT')

# SAVE IN CURRENT DIRECTORY

np.save('gdgt_sensor.npy',gdgt_proxy)

# 2.2 Carbonate Sensor

carb_proxy = carb.carb_sensor(LST,d18Ow,isoflag=-1,model='ONeil')

# SAVE IN CURRENT DIRECTORY

np.save('carb_sensor.npy',carb_proxy)

# 2.3 Leaf Wax Sensor

sample_input = 'IsoGSM_dDP_1953_2012.txt'
dDp  = np.loadtxt(sample_input) 
fC_3 = 0.7 #fraction of C3 plants
fC_4 =  0.3  #fraction of C4 plants
eps_c3 = -112.8 #pm 34.7
eps_c4 = -124.5 #pm 28.2

wax_proxy = leafwax.wax_sensor(dDp,fC_3,fC_4,eps_c3,eps_c4)

# define the error range on the epsilon (apparent fractionation) measurement:
eps_c3_err=34.7
eps_c4_err=28.2

# add uncertainties in apparent fractionation via monte-carlo resampling process:
delta_d_wax_mc,Q1,Q2=leafwax.wax_uncertainty(dDp,fC_3,fC_4,eps_c3,eps_c4,eps_c3_err,eps_c4_err)

# where Q1 is the 2.5th percentile, Q2 is the 97.5th percentile of the 1000 MC realizations

# SAVE IN CURRENT DIRECTORY

np.save('leafwax_sensor.npy',wax_proxy)

np.save('leafwax_sensor_1000mc.npy',delta_d_wax_mc)
np.save('leafwax_sensor_975CI.npy',Q2)
np.save('leafwax_sensor_25CI.npy,Q1')

print('Sensor model finished....')

#==============================================================
# L4. RUN ARCHIVE MODEL
#==============================================================

import lake_archive_compact as compact

print 'Running archive model...'

# 4.1 SEDIMENTATION

Sbar = 100.0         # Sedimentation rate, cm/kyr
S = Sbar/1000.        # Sedimentation rate, cm/yr

# Without compaction:

T = 21000.            #years !set from the envriomental model! 
H = S*T               #meters

# TODO: NOTE THAT IF YOU HAVE DEPTH MEASUREMENTS THIS SHOULD BE SPECIFIED HERE

depth_profile = H     # total depth of core [m]

h=np.linspace(0,H,num=100)
depth_profile=H # total depth of core


# Trends in porosity ~ reflect inhomogeneities in relatve compaction:
# Apply correction based on conservation of dry sediment volume wherein the thickness of each layer
# is adjusted so as to remove trends in porosity.

# Run porosity function to calculate compaction due to porosity in
# lake sediments:

phi,z = compact.porosity(depth_profile)
phi_0 = 0.95          # porosity of sediments at site (typical value for lake seds = 0.99)

np.save('phi.npy', phi)
np.save('z.npy',z)
# Now adjust the compacted sediment depth scale based on porosity

h_prime = np.zeros(len(phi)) # depth and age at each index. 
for i in range(len(phi)):
    h_prime[i] = h[i]*(1-phi_0)/(1-phi[i])

# save depth-age profile as independent output
#==============================================================

# 4.2 BIOTURBATION ~ TURBO2 (Trauth, 2013)

import lake_archive_bioturb as bio

years = np.arange(1979,2015,1) # READ THE YEARS.
age   = years
mxl   = np.ones(36)*4.
abu   = np.ones(36)*200.  # abundance of species carrier 1
abu2 = abu[-37:-1]
iso   = pseudoproxy
lngth = len(pseudoproxy)
numb  = 10

oriabu,bioabu,oriiso,bioiso = bio.bioturbation(abu2,iso,mxl,numb)

#==========================================================================
# L5. APPLY OBSERVATION MODEL
#==========================================================================
print 'Running observation model...'

import csv
from numpy import genfromtxt 
from rpy2.robjects import FloatVector
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
# require rpy2 installed

import lake_obs_bchron2 as bchron

data = genfromtxt('TEX86_cal.csv',
                delimiter = ',',
                names=True,
                dtype=None
                #mask module may be used when better understood
                #,usemask=True
                )

# example input file: user have to provide their own input

year=data['AGE']
depth=data['DP']
sds=data['SD']

r = robjects.r #robjective
calCurves=np.repeat('normal', 26)

nyears=21000. # that should fixed from the input
d=500.   # same as H in the sendimentation.

ages = FloatVector(year)    # age estimate
sd  = FloatVector(sds)# SD of ages
positions = FloatVector(depth)     # position in core in cm
calCurves = StrVector(calCurves)
predictPositions = r.seq(0,d,by=d/nyears)   # note d= total depth of core (60cm)
# Specify extractDate (top-most age of the core, in years BP. 1950=0BP)
extractDate= -55

# Call BCHRON observation model to return CHRONS (in years BP)

# NOTE: long run time waming
chronBP, depth_horizons = bchron.chron(ages,sd,positions,calCurves,predictPositions,extractDate)

# recast in years CE
chronCE = np.fliplr(1950 - chronBP)

#==========================================================================

# 3.2: Analytical Uncertainty Model:

print 'Adding uncertainty...'
#Enter final simulated speleothem record (choose from above options)
X=pseudoproxy
#X=d18O_wm2

#3.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
lake_upper, lake_lower = analytical_err_simple(X,sigma)

#3.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
lake_Xn=analytical_error(X,sigma)

#====================================================================
# Save whatever needs to be saved
print 'Saving data...'
outdir='./results/'
np.save(outdir+"lake_Xn.npy",lake_Xn)
#Whatever else...
#====================================================================

# ~
# ~
# ~
# ~

#======================================================================
# L7 PLOTTING EXAMPLES
#======================================================================

# Plotting!  Uncomment to plot compaction profile

# plt.style.use('ggplot')

# fig2 = plt.figure(1)

# fig2.add_subplot(1,2,1)

# plt.plot(z,phi,'b')
# plt.xlabel(r'Depth (m)',fontsize=16)
# plt.ylabel(r'Porosity Profile ($\phi$) (unitless)',fontsize=16)
# plt.title(r'Porosity ($\phi$) Profile in Sediment Core',fontsize=16)
# plt.grid('on',color='grey')
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.xticks(fontname = "Helvetica")  # This argument will change the font.
# plt.yticks(fontname = "Helvetica")  # This argument will change the font.

# #======================================================================

# fig2.add_subplot(1,2,2)
# plt.plot(z,h_prime,'r',label=r'Compacted Layer')
# plt.plot(z,h,'b',label=r'Non-Compacted Original Layer')
# plt.xlabel(r'Depth (m)',fontsize=16)
# plt.ylabel(r'Sediment Height (m)',fontsize=16)
# plt.title(r'Depth Scale w/Compaction in Sediment Core',fontsize=16)
# plt.grid('on',color='grey')
# plt.tick_params(axis='both', which='major', labelsize=14)
# plt.xticks(fontname = "Helvetica")  # This argument will change the font.
# plt.yticks(fontname = "Helvetica")  # This argument will change the font.
# plt.legend(loc=2)
# plt.show()


# # ======================================================================
