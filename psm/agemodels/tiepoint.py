# Sylvia Dee <sdee@usc.edu>
# PRYSM
# OBSERVATION SUB-MODEL: Tie-Point Age Models
# Function 'CHRON'
# Modified [sdee] 11/11/2014

#====================================================================

'''
    DOCSTRING: FUNCTION 'CHRON'
    Observation Model/Age Model for Tie-Point Chronologies

    Use: generate an ensemble of age-perturbed data given a series of dates.

    This function uses BCHRON, a published R-package, modified for Python after Haslett & Parnell (2008)
    See documentation for BCHRON at: <http://cran.r-project.org/web/packages/Bchron/Bchron.pdf>
    Reference: <sHaslett, J., and A. Parnell, A simple monotone process with application to
    radiocarbon-dated depth chronologies, Journal of the Royal Statistical Society: Series
    C (Applied Statistics), 57(4), 399-418, 2008.>

    INPUT ARGUMENTS:

        ages        A vector of ages (most likely 14C)
        ageSds      A vector of 1-sigma values for the ages given above
        positions   Position values (e.g. depths) for each age
        positionThicknesses (optional) Thickness values for each of the positions. By default set to zero
        calCurves   A vector of values containing either intcal13, shcal13, marine13, or normal.
                    Should be the same length the number of ages supplied. Non-standard
                    calibration curves can be used provided they are supplied in the same format
                    as those previously mentioned and are placed in the same directory. Normal
                    indicates a normally-distributed (non-14C) age.
        ids         (optional) ID names for each age
        outlierProbs        (optional) A vector of prior outlier probabilities, one for each age. Defaults to
                            0.01
        predictPositions    (optional) A vector of positions (e.g. depths) at which predicted age values are
                            required. Defaults to a sequence of length 100 from the top position to the
                            bottom position
        pathToCalCurves(optional) File path to where the calibration curves are located. Defaults to the
                    system directory where the 3 standard calibration curves are stored.
        iterations  (optional) The number of iterations to run the procedure for
                    burn (optional) The number of starting iterations to discard
        thin        (optional) The step size for every iteration to keep beyond the burnin
        extractDate (optional) The top age of the core. Used for extrapolation purposes so that no
                    extrapolated ages go beyond the top age of the core. Defaults to the current year
        maxExtrap   (optional) The maximum number of extrapolations to perform before giving up.
                    Useful for when large amounts of extrapolation are required, i.e. some of the
                    predictPositions are a long way from the dated positions
        thetaMhSd   (optional) The Metropolis-Hastings standard deviation for the age parameters
        muMhSd      (optional) The Metropolis-Hastings standard deviation for the Compound PoissonGamma mean
        psiMhSd     (optional) The Metropolis-Hastings standard deviation for the Compound PoissonGamma scale
        ageScaleVal (optional) A scale value for the ages. Bchronology works best when the ages are
                    scaled to be approximately between 0 and 100. The default value is thus 1000
                    for ages given in years.
        positionScaleVal    (optional) A scale value for the positions. Bchronology works best when the
                            positions are scaled to be approximately between 0 and 100. The default value
                            is thus 100 for positions given in cm.
    Details:
        The Bchronology function fits a compound Poisson-Gamma distribution to the increments between
        the dated levels. This involves a stochastic linear interpolation step where the age gaps are Gamma
        distributed, and the position gaps are Exponential. Radiocarbon and non-radiocarbon dates (including
        outliers) are updated within the function also by MCMC.

    KEYWORD-ARGUMENTS:
        extractDate: enter your extraction date for the core.
        (DEFAULT=-64, which is equivalent to 2014. This is the date of the core
        extraction, with years BP 1950=0).

    OUTPUT:
        chrons: ensemble of modified chronologies (1000 = default)
        depth_horizons: depths in core corresponding to each date horizon.
        (Included to facilitate plotting record with depth).

    USAGE:
        Bchronology(ages, ageSds, positions, positionThicknesses = rep(0, length(ages)),
        calCurves = rep("intcal13", length(ages)), ids = NULL,
        outlierProbs = rep(0.01, length(ages)), predictPositions = seq(min(positions),
        max(positions), length = 100), pathToCalCurves = system.file('data',
        package = "Bchron"), iterations = 10000, burn = 2000, thin = 8,
        extractDate = 1950 - as.numeric(format(Sys.time(), "%Y")), maxExtrap = 500,
        thetaMhSd = 0.5, muMhSd = 0.1, psiMhSd = 0.1, ageScaleVal = 1000,
        positionScaleVal = 100)


    Example:

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

def chron(ages,sd,positions,calCurves,predictPositions,extractDate=-64):

#======================================================================================
# Initialization
#======================================================================================

    import numpy as np
    import matplotlib.pyplot as plt
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    r = robjects.r

# Install Necessary R packages to do the analysis

    utils = importr("utils")
    utils.chooseCRANmirror(ind=1)
    packnames = ('Bchron', 'stats','graphics')
    from rpy2.robjects.vectors import StrVector
    utils.install_packages(StrVector(packnames))
    Bchron=importr('Bchron')

#======================================================================================
# Run Bchrology and store simulated chronologies as variable 'ages'
#======================================================================================

    ages = Bchron.Bchronology(ages=ages, ageSds=sd,
    positions=positions,
    calCurves=calCurves, predictPositions=predictPositions, extractDate=extractDate)

    # Plot the Ages and the 5,95% CI as a check (optional)

    r.plot(ages,main="Carlsbad Cave Age Uncertainties 1000-2005",xlab='Age (cal years BP)',ylab='Depth (mm)',las=1)

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

    return chrons, depth_horizons

