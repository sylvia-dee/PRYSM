# Sylvia Dee
# PSM d18O Cellulose
# function 'cell_sensor'
# Modified 07/02/2014 <sdee@usc.edu>

def cellulose_sensor(t,T,P,RH,d18Os,d18Op,d18Ov,flag=1.0,iso=True):

#==============================================================
    """
    DOCSTRING: Function 'cellulose_sensor'
    DESCRIPTION: Sensor Model in psm_cellulose to compute d18O of cellulose.

    ENVIRONMENTAL/GCM INPUTS:

    t       time (years, vector array)
    T       (Temperature, K)
    P       (Precipitation, mm/day)
    RH      (Relative Humidity, %)
    d18Os   (Isotope Ratio of Soil Water)
    d18Op   (Isotopes Ratio of Precipitation)
    d18Ov   (Isotope Ratio of Ambient Vapor at Surface Layer)

    OUTPUT:

    d18O(Cellulose) [permil]

    WORKING VARIABLES:

    dlw       = d18O leaf water
    d18longts = d18O precip/soil water
    dwvlongts = d18O vapor
    droden    = d18O tree cellulose

    Notes: d18O values of leaf cellulose are 27 elevated above leaf water (Sternberg, 1989; Yakir and DeNiro, 1990)~ = ek
    The mechanism for the 27 enrichment is the carbonyl\u2013water interaction during biosynthesis (Sternberg and DeNiro, 1983).
    Because the fractionation factor is the same for autotrophic and heterotrophic fractionation of oxygen,
    there is no need to distinguish between the two.

    Two models are available: Roden, Evans

    1. Roden (Roden et al. 2003): (NOTE: FLAG = 0 for Roden)
    2. Evans (Evans et al. 2007): (NOTE: FLAG = 1 for Evans)
    """
#==============================================================
# 1. Roden Model
#==============================================================
    import numpy as np

    if flag == 0:
# 1.1 Define Constants
        ek = -27.
        ecell = 27.
        ee = 9.8
# 1.2 Load Relative Humidity
        relhum = RH/100.
        fo = 0.35
# 1.3 Compute delta value of Leaf Water [dlw] (Ciasis, 2000)
        d18longts = d18Os
        dwvlongts = d18Ov
        dlw=ee + (1.0-relhum)*(d18longts-ek)+relhum*(dwvlongts)
# 1.4 Compute d18O of Tree Cellulose ~ Yearly averages (timeseries) for 100 years

        droden=fo*(d18longts+ecell) + (1.0-fo)*(dlw+ecell)
        dcell=droden
        return dcell

#==============================================================
# 2. Evans Model
#==============================================================

    elif flag == 1.0:
        # Convert Units for Inputs (USER CHECK)
        ta = T-273.15
        pr = P*30
        rh = RH
        # Parameters for Evans Model (see Evans 2007, doi:10.1029/2006GC001406)
        dfs=32.
        dfb=21.
        efc=27.
        blc=2100.
        pl=0.018
        feo=0.42
        fwxm= 1.0
        c1 = 0.74
        c2 = -0.0285
        c3 = 10.
        c4 = 23./33.
        c5 = 3.
        c6 = 6.
        c7 = 27.17
        c8 = 0.45
        c9 = 20.
        c10= 4290.
        c11=-43.
        #==============================================================
        # Sensor Model
        #==============================================================
        # Option 1: Use Isotope-Enabled Model Output for dS, dV
        if iso==True:
            o18sw = d18Os
            lt= c3+c4*ta       # leaf temperature
            tc=-c6*c5+273.15+ta
            alv=np.exp(1137.0/tc**2.0-0.4156/tc-0.00207)
            o18wva = d18Ov
        # Option 2: use Precipitation to calculate soil water, parameterize vapor:
        else:
            o18sw = c1 + c2*pr    # soil water
            lt= c3+c4*ta       # leaf temperature
            tc=-c6*c5+273.15+ta
            alv=np.exp(1137.0/tc**2.0-0.4156/tc-0.00207)
            o18wva=1000*(1./alv-1)+o18sw

            # Majoube, 1971 referenced in Gonfiantini et al., 2001:
            # equilibrium fraction between water liquid and vapor
            # o18wva = o18sw - 8

        # Set Constants

        rsmow=0.0020052     # 18/16 ratio for intl working std
        cw=55500.           # concentration of water (mol/m^3)
        dw=0.00000000266    # diffusivity of H218O in water (m^2/s)

        #Environmental Conditions

        svpa=(6.13753*np.exp(ta*((18.564-(ta/254.4)))/(ta+255.57))) # saturated vapour pressure of the air (mbar)
        avpa=(rh/100.)*svpa             # actual vapour pressure of the air (mbar)
        rsw=((o18sw/1000.)+1.)*rsmow    # 18O/16O of source water
        rwv=((o18wva/1000.)+1.)*rsmow   # 18O/16O of water vapour
        do18wvs=((rwv/rsw)-1.)*1000.    # delta 18O of water vapour (per mil, wrt source water)

        import scipy as SP

        # leaf evaporative conditions
        vpl=(6.13753*SP.exp(lt*((18.564-(lt/254.4)))/(lt+255.57))) # leaf internal vapour pressure (mbar)
        vpd=vpl-avpa # leaf - air saturation vapor pressure deficit (mbar)
        sc=c10/vpd + c11 # Lohammar function for sc=f(vpd) vpd units in mb sc units in mmol/m^2/sec
        tr=sc*vpd/1013.0 # tr = sc*vpd/Patm (mb)
        # Note for tr calculation, atmospheric pressure is assumed constant at sea level value
        # this should be adjusted for sites at appreciable elevation, or input as monthly timeseries
        pe=((tr/1000.)*pl)/(cw*dw) # Peclet number
        # fractionation factors
        tdf=(((1./sc)*dfs)+((1./blc)*dfb))/((1./blc)+(1./sc)) # total diffusive fractionation (per mil)

        a1 = (lt+273.16)**2
        a2=(0.4156/(lt+273.16))
        a3=0.0020667

        evpf=(SP.exp((1137./a1)-a2-a3)-1)*1000; # equilibrium vapour pressure fractionation (per mil)

        do18e=(((1.+(evpf/1000.))*(1.+(tdf/1000.)+(((do18wvs/1000.)-(tdf/1000.))*(avpa/vpl))))-1.)*1000. # Craig-Gordon evaporation-site water (D\18Oe) (per mil wrt source water )
        do18l=(do18e*(1.-(SP.exp(-pe))))/pe # Bulk leaf water (D18OL) (per mil wrt source water)
        do18s=do18l+efc # Sucrose D18O (per mil wrt source water)
        do18c=do18l*(1.-(fwxm*feo))+efc # Cellulose D18O (per mil wrt source water)

        #convert back into ratios wrt SMOW:
        M=1000.
        o18e=(((do18e/M+1.)*rsw)/rsmow-1.)*M
        o18l=(((do18l/M+1.)*rsw)/rsmow-1.)*M
        o18s=(((do18s/M+1.)*rsw)/rsmow-1.)*M
        o18c=(((do18c/M+1.)*rsw)/rsmow-1.)*M

        # correct to observed o18 given observed mean (c7)
        # and fraction (c8) water from meteoric sources
        o18cc=o18c
        o18cc=c8*(o18c-(np.mean(o18c)-c7))+(1-c8)*c7
        dcell=o18cc

        return dcell
