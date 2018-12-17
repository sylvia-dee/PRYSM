#!/usr/bin/env python
# Sylvia Dee
# PSM d18O Ice Cores
# ENVIRONMENT SUB-MODEL
# a.k.a. driver script
# Modified 01/27/2015 <sdee@usc.edu>

"""
    Ice Core d18O Environment Model (Wrapper Script)
    ~ Load Environmental Variables to compute d18O calcite/dripwater.

    INPUTS:
        Isotope-Enabled GCM INPUTS:
        T       (Temperature, K)
        P       (Pressure, Atm)
        d18O    (Precipitation, permil)

    Wrapper script runs the following submodels:

    icecore_sensor  (weighted d18O precip, altitude and temperature bias corrections)
    icecore_archive (diffusion, compaction)
    icecore_bam     (time uncertainties, analytical uncertainties)

    OUTPUTS:
        d18O Ice Core Timeseries
"""
#======================================================================
# E0. Initialization
#======================================================================
import numpy as np
import scipy
from scipy import integrate
import matplotlib.pyplot as plt
import math
from math import exp
from math import sqrt
from math import *

from psm.icecore.sensor import icecore_sensor
from psm.icecore.archive import icecore_diffuse
from psm.agemodels.banded import bam_simul_perturb
from psm.aux_functions.analytical_error import analytical_error
from psm.aux_functions.analytical_err_simple import analytical_err_simple

#======================================================================
# E1. Load input data/all environmental variables (annual averages)
#======================================================================

# LOAD TESTER FILES (OPTIONAL) or input your own data here.
datadir='test_data_icecore/'
print('Loading data from ',datadir)
# Near-Surface Air Temp [K]
T=np.load(datadir+'temperature.npy')

# Accumulation Rate [m]
accum=np.load(datadir+'accumulation.npy')

# Core Depth
depth=np.load(datadir+'total_icecore_depth.npy')

# Depth horizons (accumulation per year corresponding to depth moving down-core)
depth_horizons=np.load(datadir+'depth_horizons.npy')

#delta-18-O of precipitation (weighted) [permil]
d18Op=np.load(datadir+'d18O_P.npy')

# Mean Surface Pressure [Atm]
P=1.0

# Set Time Axis
time=np.arange(1000,2005,1)

#======================================================================
# E2. Calculate Precip-Weighted Delta Values at one location
#======================================================================

deltaP = d18Op#your precip isotope data would go here

# reshape data sets: singular vector
# example: deltaP=deltaP.reshape(1005)

#======================================================================
# E3. CALL SENSOR MODEL
#======================================================================
print('Running sensor model...')
# 4. Apply icecore_sensor to extract precipitation-weighted d18o record for each core
#    and compute altitude, temperature corrections. (Please see docstring icecore_sensor).

d18O     = deltaP      # your dataset loaded here
alt_diff = 3524.0      # alt_diff at location (m) (THIS IS FOR SPEEDY-QUELCCAYA)
d18Oice  = icecore_sensor(time,d18O,alt_diff)

# returns: icecore

#======================================================================
# E4. CALL ARCHIVE MODEL
#======================================================================
print('Running archive model...')
# This archive model will calculate diffusion and compaction
# (Please see docstrings: diffusivity, icecore_diffuse)

# NOTE: tester file has accumulation in meters per year. Below is optional unit conversion.

#accum=accum*365.0       # multiple by 365 days to get yearly accumulation in mm
#b = accum/1000.         # convert mm/yr to m/yr, accumulation rate (e.g. 1.3 m/year)

b=accum
core_length=np.cumsum(b)
depth = core_length[-1]
#depth_horizons=np.diff(core_length)#enter depth horizons downcore

# Set inputs for Archive Model

# For the Density Profile, we use average T, b
T=np.mean(T)
b=np.mean(b)

# Optional: set manually
#P = 1.# pressure in atmospheres
#T = 296.02# temperature (K), default value = 260 K

dz = min(depth_horizons)/10.#or 0.06# stepping in depth (m)
drho = 0.5# step in density (kg/m^3)
depth= depth# enter your ice core total depth/length in meters

# Run the archive model:
z, sigma, D, time_d, diffs, ice_diffused, rho, zieq = icecore_diffuse(d18Oice,b,time,T,P,depth,depth_horizons, dz,drho)

#======================================================================
# E5. CALL OBSERVATION MODEL
#======================================================================
print('Running observation model...')
# 5.1 Specify and model rate of annual layer miscount

X = ice_diffused
X = X.reshape(len(X),1)
t = time
tp, Xp, tmc=bam_simul_perturb(X,t,param=[0.01,0.01],name='poisson',ns=1000,resize=0)
#======================================================================

# 5.2: Analytical Uncertainty Model:
print('Adding uncertainty...')
#5.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
ice_upper, ice_lower = analytical_err_simple(X,sigma)

#5.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
ice_Xn=analytical_error(X,sigma)

#====================================================================
# Save whatever needs to be saved
print('Saving data...')
outdir='./results/'
np.save(outdir+"ice_Xn.npy",ice_Xn)
np.save(outdir+"ice_diffused.npy",ice_diffused)
np.save(outdir+"ice_depth.npy",z)
np.save(outdir+"ice_time_d.npy",time_d)
#Whatever else...
#====================================================================
