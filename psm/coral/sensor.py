# Sylvia Dee <sdee@usc.edu>
# PSM d18O Coral Aragonite
# SENSOR MODEL
# Function 'pseudocoral'
# Modified 10_16_2015 <sdee@usc.edu>
# Modified 11_17_2015 <sylvia_dee@brown.edu>
#====================================================================

def pseudocoral(lat,lon,SST,SSS,d18O=-1,species="default",b1=0.3007062,b2=0.2619054,b3=0.436509,b4=0.1552032,b5=0.15):

    """
        DOCSTRING: Function 'pseudocoral' produces a d18O-coral record given SST, SSS, and global position.
           The function is based on the forward model published by [Thompson, 2011]:
           <Thompson, D. M., T. R. Ault, M. N. Evans, J. E. Cole, and J. Emile-Geay (2011),
           Comparison of observed and simulated tropical climate trends using a forward
           model of coral \u03b418O, Geophys.Res.Lett., 38, L14706, doi:10.1029/2011GL048224.>
           Returns a numpy array that is the same size and shape as the input vectors for SST, SSS.
        
        Input parameters:
            Latitude    [lat]       [-90, 90]
            Longitude   [lon]       [0, 360]
            SST         [SST]       SST ANOMALY Units = degrees Celsius
            SSS         [SSS]       SSS ANOMALY Units = PSU
            Please note that input files must be read in as 1-D vectors.
            
        Output:
            Returns a Numpy Array, ['coral'], which is saved in the main program.
            Optional Inputs: Please note that the function will use default params for
            d18O vs. SSS, a, and b unless specified and placed in the appropriate location as shown below.
            delta_18O [permil]: If the user does not have an input field for d18O, you must put -1 in its 
            position in the call and use the equation for SSS below.

            **Note: [always use d18O (seawater) if available!]**

            The user may specify values for parameters a and b in the coral d18O forward model,
            which are currently set to the following default values:

            a = coral - SST regression, accounting for species-dependent variations in the SST-18O dependence.
            The code assigns this value universally at -.22, but published values vary by species/genus (i.e. Moses, 2006)
            Specify keyword input 'species=' to set slope:
            species = 
                    "Default":      a = -0.22
                    "Porites_sp":   a = -.26178
                    "Porites_lob":  a = -.19646
                    "Porites_lut":  a = -.17391
                    "Porites_aus":  a = -.21
                    "Montast":      a = -.22124
                    "Diploas":      a = -.14992

            b = d18O-SSS slope, as defined by location in: Legrande & Schmidt, (2006)
                    <LeGrande, A. N., and G. A. Schmidt (2006), Global gridded data set of the oxygen isotopic          composition in seawater, Geophys. Res. Lett., 33, L12604, doi:10.1029/2006GL026011.>
                    b1 = 0.31   [Red Sea]
                    b2 = 0.27   [Tropical Pacific]
                    b3 = 0.45   [South Pacific]
                    b4 = 0.16   [Indian Ocean]
                    b5 = 0.15   [Tropical Atlantic Ocean/Caribbean]
        Example of use: Call function 'pseudocoral' to calculate d18O_coral timeseries, lacking d18O sw
        but specifying own numeric parameters for a, b.
        pseudocoral(lat,lon,SST,SSS,-1,species="default",b1,b2,b3,b4,b5)
    """
#====================================================================



# Define slope of coral response to SST, a
    #Also, have you considered accounting for species-dependent variations in the a1 value (the SST-18O dependence)? 
    #The code assigns this value universally at -.22, but published values vary by species/genus (i.e. Moses, 2006). 
    #I've copied and pasted a chunk of code from one of my scripts that does this brute-force... Hope it helps.

    if species == "Porites_sp":
        a=-.26178
    elif species == "Porites_lob":
        a=-.19646
    elif species == "Porites_lut":
        a=-.17391
    elif species == "Porites_aus":
        a=-.21
    elif species== "Montast":
        a=-.22124
    elif species == "Diploas":
        a=-.14992
    elif species == "default":
        a=-0.22
    else:
        a=-0.22

# Define SSS-d18O(sw) slopes by region. Note that we multiply all slopes by VPDB/VSMOW [0.0020052/0.0020672] ~ 0.97 to correct d18O values.

# Because d18Osw is reported relative to VSMOW, we have to convert the slope to VPDB so
# that the units of the forward model output are in VPDB (and thus directly comparable
# to published d18O records, in units of VPDB). The way this affects the forward model is
# in the d18Osw-SSS slope.

    V_corr = 0.97002

    slope1 = V_corr*b1      # Red Sea
    slope2 = V_corr*b2      # Tropical Pacific
    slope3 = V_corr*b3      # South Pacific
    slope4 = V_corr*b4      # Indian Ocean
    slope5 = V_corr*b5      # Tropical Atlantic

# Given Lat/Lon pair, assign correct value to 'b'
    if lon>=32.83 and lon<=43.5 and lat>=12.38 and lat<=28.5:    # Red Sea
        b=slope1
    elif lon<=120.:              # Indian Ocean
        b=slope4
    elif lon>=270. and lat>=10.:    #Tropical Atlantic ~300 E longitude 
        b=slope5
    elif lat> -5. and lat<=13.:   # Tropical Pacific
        b=slope2
    elif lat<= -5.:              # South Pacific
        b=slope3
    else:
        b=slope2            # Default = Trop. Pacific.

#====================================================================

# Form an array of pseudocoral data

    if d18O==-1:
        coral = a*SST+b*SSS
    else:
        coral = a*SST + d18O

#====================================================================

    return coral
