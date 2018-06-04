# Sylvia Dee <sylvia@ig.utexas.edu>
# PRYSM
# PSM for Lacustrine Sediments
# OBSERVATION MODEL for age model error analysis
# Function 'lake_obs_bchron'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
# Modified 04/2/2018 <sylvia@ig.utexas.edu>

'''
    BCHRON Example:
        from rpy2.robjects import FloatVector
        from rpy2.robjects.vectors import StrVector
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        r = robjects.r
        # aply chronological errors
        d=60.                                       # total length of core in cm
        ages = FloatVector([-55,55,236,432,835])    # age estimate
        sd  = FloatVector([1,22,50,13,25])# SD of ages
        positions = FloatVector([0,5,15,35,55])     # position in core in cm
        calCurves = StrVector(['normal','normal','normal','normal','normal'])
        predictPositions = r.seq(0,d,by=d/nyears)   # note d= total depth of core (60cm)
        # Specify extractDate (top-most age of the core, in years BP. 1950=0BP)
        topDate= -55
        # Call BCHRON observation model to return CHRONS (in years BP)
        chronBP, depth_horizons = chron.Bchron(ages,sd,positions,calCurves,predictPositions,topDate)
    '''
#======================================================================================
# Install Necessary R packages to do the analysis

import numpy as np
import matplotlib.pyplot as plt
from rpy2.robjects import FloatVector
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
r = robjects.r

utils = importr("utils")
utils.chooseCRANmirror(ind=1)
packnames = ('Bchron', 'stats','graphics')
from rpy2.robjects.vectors import StrVector
utils.install_packages(StrVector(packnames))
Bchron=importr('Bchron')


#apply chronological errors
#LOAD CSV of AGES/SDs

cd LAKEPSM

import csv
from numpy import genfromtxt

data = genfromtxt('TEX86_cal.csv',
                delimiter = ',',
                names=True,
                dtype=None
                #mask module may be used when better understood
                #,usemask=True
                )
year=data['AGE']
depth=data['DP']
sds=data['SD']


calCurves=np.repeat('normal', 26)

nyears=21000.
d=500.        
ages = FloatVector(year)    # age estimate
sd  = FloatVector(sds)# SD of ages
positions = FloatVector(depth)     # position in core in cm
calCurves = StrVector(calCurves)
predictPositions = r.seq(0,d,by=d/nyears)   # note d= total depth of core (60cm)
# Specify extractDate (top-most age of the core, in years BP. 1950=0BP)
extractDate= -55
#======================================================================================
# Initialization
#======================================================================================
#======================================================================================
# Run Bchrology and store simulated chronologies as variable 'ages'
#======================================================================================

ages = Bchron.Bchronology(ages=ages, ageSds=sd,
positions=positions,
calCurves=calCurves, predictPositions=predictPositions, extractDate=extractDate)

# Plot the Ages and the 5,95% CI as a check (optional)

r.plot(ages,main="Tanganyika Age Uncertainties, LGM-Present",xlab='Age (cal years BP)',ylab='Depth (mm)',las=1)

#======================================================================================
# Extract individual chronologies
#======================================================================================

    # The variable 'theta' is the posterior estimated values of the ages, and is stored as
    # the first matrix of the List Vector class 'ages'

theta = ages[0]
theta=np.array(theta)

thetaPredict = ages[4]
thetaPredict=np.array(thetaPredict)

# Save depth horizons
depths=np.array(predictPositions)
depth_horizons=depths[:-1]

# save all the chronologies as independent time axes and plot original data with
# varying X values.

thetaPredict.shape
chrons=thetaPredict[:,:-1]      # reshape to be same size as input data
chrons.shape
#======================================================================================
# Return all chronologies
#======================================================================================
