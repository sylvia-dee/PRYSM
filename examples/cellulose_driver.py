#!/usr/bin/env python
# Sylvia Dee
# PSM d18O Cellulose
# Driver Script
# ENVIRONMENT SUB-MODEL
# Modified 01/27/15 <sdee@usc.edu>
#==============================================================
"""
Tree Cellulose ~ Load Environmental Variables to compute d18O cellulose.
This is a driver script. It is set up to walk the user through the necessary
commands to use the Cellulose PSM. The user can load there own data in the fields below,
or use the provided test data.

SPEEDY/GCM INPUTS:

    T       Temperature (K)
    d18O    (Soil Water)
    d18O    (Ambient Vapor)
    d18O    (Precipitation)
    RH      (Relative Humidity)

This driver script calls the submodels:

    cellulose_sensor
    cellulose_archive
    cellulose_obs

The final output is d18O of cellulose at a given cite, with 1000 plausible chronologies
given a specified rate of miscount.

"""
#==============================================================
# DRIVER SCRIPT
#==============================================================

# 0.0 Initialization
import numpy as np
import psm.cellulose.sensor as cs
from psm.agemodels.banded import bam_simul_perturb
from psm.aux_functions.analytical_error import analytical_error
from psm.aux_functions.analytical_err_simple import analytical_err_simple
#==============================================================
# 1.0 Load Test Variables or your own data here:
#==============================================================

# Example data directory
datadir='test_data_cellulose/'
print('Loading data from ', datadir,' ...')

# Soil Water Isotope Ratio [permil]
d18Os=np.load(datadir+'d18O_soil.npy')
d18Os=d18Os[:,1]

# Precip [mm/day]
P=np.load(datadir+'precipitation.npy')

# Near-Surface Relative Humidity [%]
RH=np.load(datadir+'relative_hum.npy')

# Ambient Vapor Isotope Ratio d18Ov [permil]
# NOTE: Take out the Lowest (surface) Level only.
d18Ov=np.load(datadir+'d18O_vapor.npy')
d18Ov=d18Ov[:,1]

# Annual Average Surface Temperature
T=np.load(datadir+'temperature.npy')

# Precipitation Isotope Ratio [permil]
d18Op=np.load(datadir+'d18O_precip.npy')

#==============================================================
# Set time Axis
#==============================================================
time = np.arange(1000,2005,1)
#==============================================================
# Check dimensions: must be 1-D Vector
#==============================================================
# example: >VAdO18_surface.shape

#==============================================================
# 2.0 Set geographic location.
#==============================================================

# Evans Model is suitable for the entirety of the tropics (30N-30S)
# Choose location : (test files are for La Selva, Costa Rica)

# Note: you can use cdutil package to pull this from a climate model field:
# example:
# >import sys, os, cdtime, cdutil, cdms2, vcs, MV2, time, datetime
# >las = cdutil.region.domain(latitude=(10.,10.),longitude = (-84.,-84.))
# >var = GLOBALVAR(las)

#==============================================================
# 3. RUN SENSOR MODEL
#==============================================================

# Call function pseudocell to compute d18O cellulose field
# (Please see docstring for cellulose_sensor for instructions).

# UNITS: T[K], P[MM/DAY], RH [%], isotope fields [permil]
# KWARGS:   [flag = 1 for Evans Model, 0 for Roden Model]
#           [iso=true to use isotope fields, false to use T, P fields]
print('Running sensor model...')
dcell = cs.cellulose_sensor(time,T,P,RH,d18Os,d18Op,d18Ov,flag=1.0,iso=True)

#==============================================================
# 4. RUN ARCHIVE MODEL
#==============================================================
# NOTE: PlACEHOLDER.

#dcell = cell_archive()

#==============================================================
# 5. RUN OBSERVATION MODEL
#==============================================================

# Call function bam_simul_perturb to compute ensemble of plausible
# age-uncertain records and account for analytical uncertainty.

# 5.1: Ring Miscount Uncertainty Model:
print('Running observation model...')
X=dcell
t=time
X=X.reshape(len(X),1)
tp, Xp, tmc = bam_simul_perturb(X,t,param=[0.05,0.05],name='poisson',ns=1000,resize=0)
#==============================================================

# 5.2: Analytical Uncertainty Model:
print('Adding uncertainty...')
#5.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
dcell_upper, dcell_lower = analytical_err_simple(X,sigma)

#5.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
dcell_Xn=analytical_error(X,sigma)
#====================================================================
# Save whatever needs to be saved
print('Saving data...')
outdir='./results/'
np.save(outdir+"dcell_Xn.npy",dcell_Xn)
#Whatever else...
#====================================================================