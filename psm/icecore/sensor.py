# Sylvia Dee
# PSM d18O Ice Cores
# SENSOR SUB-MODEL
# Modified 07/02/2014 <sdee@usc.edu>

def icecore_sensor(time,d18O,alt_diff=0.):
    '''
    DOCSTRING: FUNCTION icecore_sensor

    DESCRIPTION:
    The ice core sensor model calculates precipitation-weighted
    del18OP (i.e. isotope ratio is weighted by the amount of precipitation that accumulates)
    and corrects for temperature and altitude bias between model and site ([Yurtsever, 1975],
    0.3/100m [Vogel et al., 1975]).

    CITATION: Yurtsever, Y., Worldwide survey of stable isotopes in precipitation., Rep. Isotope Hydrology Section, IAEA, 1975.

    INPUTS:
    time: (calendar year)
    d18O: d18O of precipitation (permil)
    alt_diff: Actual Altitude-Model Altitude (meters)

    RETURNS: 'icecore'

    EXAMPLE USE:
    Apply icecore_sensor to extract precipitation-weighted d18o record for each core
    and compute altitude, temperature corrections. (Please see docstring icecore_sensor).

    d18O     = ice      # your dataset loaded here

    Enter your altitude correction in meters: (+ if too low. e.g. 500 m too low, alt_diff = 500).
    alt_diff = 500      # alt_diff at location (m)

    (Example: compute altitude correction for the core)
    z1=5670-z(QUEL)

    d18Oice = icecore_sensor(time,d18O,alt_diff)

    '''
#====================================================================
# S1. Correct for altitude bias
#====================================================================
# Altitude Effect: cooling and precipitation of heavy isotopes.
# O18~0.15 to 0.30 permil per 100m.
    alt_eff  = -0.25
    alt_corr = (alt_diff/100.)*alt_eff
    ice = d18O + alt_corr
#===================================================================
# S2. Fill single array with ice core values:
#===================================================================
    #icecore = ice.reshape(len(t))
    return ice
