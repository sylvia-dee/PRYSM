#!/usr/bin/env python
# Sylvia Dee <sdee@usc.edu>
# Modified 4/25/13
# psm.coral.env
#====================================================================

"""
CORAL PSM DRIVER SCRIPT

Data preparation script to generate simulated proxy data.
This file loads input files and environmental variables and calls
the forward model functions (sub-models).

This wrapper script prepares forward model environmental input files, which
must contain spatial information (lat, lon) and puts them into
an appropriate format that the function handles (vectors only)!

Please read docstring carefully for each function call before use.

# CORALS: calculating d18O anomalies from SST and SSS anomalies.
"""
#====================================================================
# 0. Initialization
#====================================================================
import numpy as np
import matplotlib.pyplot as plt
from psm.coral.sensor import pseudocoral
from psm.agemodels.banded import bam_simul_perturb
from psm.aux_functions.analytical_error import analytical_error
from psm.aux_functions.analytical_err_simple import analytical_err_simple

#====================================================================
# 1. LOAD CLIMATE FIELDS
#====================================================================

#Data directory
datadir='test_data_coral/' #Don't forget the /
print('Loading data from ', datadir,' ...')

# Load SST anomalies [K] (NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!)
# yearly
ssta=np.load(datadir+'Palmyra_SST_Anomalies_Yearly_850-1850.npy')
# monthly
ssta_m=np.load(datadir+'Palmyra_SST_Anomalies_Monthly_850-1850.npy')

# Load SSS anomalies [psu] (NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!)
# yearly
sssa=np.load(datadir+'Palmyra_SSS_Anomalies_Yearly_850-1850.npy')
# monthly
sssa_m=np.load(datadir+'Palmyra_SSS_Anomalies_Monthly_850-1850.npy')

# set time axis

time=np.arange(850,1850,1)
#====================================================================
# 2. DATA PREP: Exception Handling
#====================================================================
print('Preparing data...')
# 2.1 Convert Lats and Lons to standard format (if they aren't already).

# lon: convert (-180: +180) to (0: 360)
# lat: convert (0: 180) to (-90: +90)

# Enter your coordinates here:
lon=197.92
lat=5.8833

# make sure there are no negative-longitude coordinates: [0 to 360] only.
# longitude
if (lon<0.):
    lon = lon+360.

# latitude
if (lat>90.):
    lat = lat-90.

# 2.2 Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
sst=ssta #AAA
sss=sssa #AAA
temp_flag = any(sst>200)

for i in range(len(sst)):
    if (temp_flag):
        sst[i] = sst[i]-274.15

#====================================================================
# 3. CALL SENSOR MODEL: Call function 'pseudocoral' (see doc) to compute pseudocoral
#    timeseries and spatial field.
#====================================================================

# NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!
print('Running sensor model...')
coral = np.zeros(len(time))  # this will initialize a [Time x Lat x Lon] matrix of coral values.

# Fill coral array with data same size as input vectors.
# Important: if you have d18O of seawater, add a 5th argument 'd18O' and input data as vector array (default is -1).
# Inputs: lat, lon, SST, SSS OR d18O
# Optional d18O/T slope parameters:
# a= -0.22,b1=0.3007062,b2=0.2619054,b3=0.436509,b4=0.1552032 (see doctring)

for i in range(len(time)):
    coral[i] = pseudocoral(lat,lon,sst[i],sss[i])

#======================================================================
# 4. CALL OBSERVATION MODEL
#======================================================================

# 4.1 Specify and model rate of annual layer miscount: BAM (see doctring)
print('Running observation model...')
X = coral
X = X.reshape(len(X),1)
tp, Xp, tmc=bam_simul_perturb(X,time,param=[0.02,0.02],name='poisson',ns=1000,resize=0)
#======================================================================
# 4.2: Analytical Uncertainty Model:

#4.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
coral_upper, coral_lower = analytical_err_simple(X,sigma)

#4.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
coral_Xn=analytical_error(X,sigma)
#====================================================================
# Save coral timeseries fields as numpy arrays in current directory.
print('Saving time series...')
outdir='./results/'
np.save(outdir+"simulated_coral_d18O.npy",coral)
np.save(outdir+"coral_age_perturbed.npy",Xp)
#coral_error_bounds=np.save('coral_error.npy',coral_upper, coral_lower)
#====================================================================
