# Sylvia Dee
# PSM d18O Cellulose
# function 'cell_sensor'
# Modified 07/02/2014 <sdee@usc.edu>

import numpy as np
import scipy as SP
import six


class CelluloseSensor(object):
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

    1. Roden (Roden et al. 2003): (NOTE: model_flag = 0 or = 'roden' for Roden)
    2. Evans (Evans et al. 2007): (NOTE: model_flag = 1 or = 'evans' for Evans)
    """
    def __init__(self, t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='roden', iso=True):

        self.t = t
        self._T = T
        self._P = P
        self._RH = RH
        self._d18Os = d18Os
        self._d18Op = d18Op
        self._d18Ov = d18Ov
        self._model_flag = model_flag
        self._iso = iso

        # sanitation/sanity checks checks go here
        if self.t.ndim > 2:
            ValueError('time must be a two dimensional array')
        if self._T < 273:
            Warning('bad temperature given, less than absolute zero')
        if self._P < 0:
            Warning('bad precipiation data given, less than zero')
        if isinstance(self._model_flag, six.string_types):
            self._model_flag = self._model_flag.lower() # make it all lowercase for convenience
        # etc....

        # initialize the appropriate model
        if self._model_flag in [0, 'roden']:
            self._model = self._RodenModel(self._RH, self._d18Os, self._d18Ov)
        elif self._model_flag in [1, 'evans']:    
            self._model = self._EvansModel(self._T, self._P, self._RH, self._d18Os,
                                          self._d18Ov, self._iso)
        else:
            raise ValueError('invalid model selection')


    def run_sensor(self):
        """
        After the sensor has been initialized, then we can run it
        """
        self._model.run()
        self.dcell = self._model.dcell


    def set_RH(self, value): 
        """ 
        Allows user to change relative humidity to the recompute model 
        without needing to reinitialize the entire model.
        """ 
        self._model.RH = value


    def set_d18Os(self, value): 
        """ 
        Allows user to change soil moisture isotope ratio to recompute model 
        without needing to reinitialize the entire model. 
        """ 
        self._model.d180s = value


    def set_d18Ov(self, value): 
        """ 
        Allows user to change vapor isotope ratio to recompute model 
        without needing to reinitialize the entire model. 
        """ 
        self._model.d180s = value


    class _RodenModel(object):
        """
        Roden model from (Roden et al. 2003)
        """
        def __init__(self, RH, d18Os, d18Ov):
            # 1.1 Define Constants
            self.ek = -27.
            self.ecell = 27.
            self.ee = 9.8
            
            # 1.2 Load Relative Humidity
            self.RH = RH / 100.
            self.d18Os = d18Os
            self.d18Ov = d18Ov
            self.fo = 0.35


        def run(self):
            # 1.3 Compute delta value of Leaf Water [dlw] (Ciasis, 2000)
            self.d18longts = self.d18Os
            self.dwvlongts = self.d18Ov

            dlw = self.ee + (1.0-self.RH) * (self.d18longts-self.ek) + self.RH * (self.dwvlongts)

            # 1.4 Compute d18O of Tree Cellulose ~ Yearly averages (timeseries) for 100 years
            droden = self.fo * (self.d18longts+self.ecell) + (1.0-self.fo)*(dlw+self.ecell)
            self.dcell = droden

            # we also return dcell as a convenience if desired
            return self.dcell
        


    #==============================================================
    # 2. Evans Model
    #==============================================================

    class _EvansModel(object):
        """
        Evans model from (Evans et al. 2007)
        """
        def __init__(self, T, P, RH, d18Os, d18Ov, iso):
            # Convert Units for Inputs (USER CHECK)
            self.ta = T - 273.15
            self.pr = P * 30
            self.RH = RH
            self.d18Os = d18Os
            self.d18Ov = d18Ov


        def run(self):

            # Parameters for Evans Model (see Evans 2007, doi:10.1029/2006GC001406)
            dfs = 32.
            dfb = 21.
            efc = 27.
            blc = 2100.
            pl = 0.018
            feo = 0.42
            fwxm = 1.0
            c1 = 0.74
            c2 = -0.0285
            c3 = 10.
            c4 = 23. / 33.
            c5 = 3.
            c6 = 6.
            c7 = 27.17
            c8 = 0.45
            c9 = 20.
            c10 = 4290.
            c11 = -43.

            lt = c3 + c4 * self.ta # leaf temperature
            tc = -c6 * c5 + 273.15 + self.ta
            alv = np.exp(1137.0 / tc**2.0 - 0.4156 / tc - 0.00207)            
            if iso:
                # Option 1: Use Isotope-Enabled Model Output for dS, dV
                o18sw = self.d18Os
                o18wva = self.d18Ov
            
            else:
                # Option 2: use Precipitation to calculate soil water, parameterize vapor:
                o18sw = c1 + c2 * self.pr    # soil water
                o18wva=1000*(1./alv-1)+o18sw

                # Majoube, 1971 referenced in Gonfiantini et al., 2001:
                # equilibrium fraction between water liquid and vapor
                # o18wva = o18sw - 8

            # Set Constants
            rsmow = 0.0020052     # 18/16 ratio for intl working std
            cw = 55500.           # concentration of water (mol/m^3)
            dw = 0.00000000266    # diffusivity of H218O in water (m^2/s)

            # Environmental Conditions
            svpa = (6.13753*np.exp(self.ta*((18.564-(self.ta/254.4)))/(self.ta+255.57))) # saturated vapour pressure of the air (mbar)
            avpa = (rh/100.)*svpa             # actual vapour pressure of the air (mbar)
            rsw = ((o18sw/1000.)+1.)*rsmow    # 18O/16O of source water
            rwv = ((o18wva/1000.)+1.)*rsmow   # 18O/16O of water vapour
            do18wvs = ((rwv/rsw)-1.)*1000.    # delta 18O of water vapour (per mil, wrt source water)


            """
            Is there a reason for the switch below to using scipy's exp 
            instead of numpys exp, as you did above? Doesn't matter, 
            but no need to load scipy if not needed...
            """

            # Leaf evaporative conditions
            vpl = (6.13753*SP.exp(lt*((18.564-(lt/254.4)))/(lt+255.57))) # leaf internal vapour pressure (mbar)
            vpd = vpl-avpa # leaf - air saturation vapor pressure deficit (mbar)
            sc = c10/vpd + c11 # Lohammar function for sc=f(vpd) vpd units in mb sc units in mmol/m^2/sec
            tr = sc*vpd/1013.0 # tr = sc*vpd/Patm (mb)
            # Note for tr calculation, atmospheric pressure is assumed constant at sea level value
            # this should be adjusted for sites at appreciable elevation, or input as monthly timeseries
            pe = ((tr/1000.)*pl)/(cw*dw) # Peclet number
            # fractionation factors
            tdf = (((1./sc)*dfs)+((1./blc)*dfb))/((1./blc)+(1./sc)) # total diffusive fractionation (per mil)

            a1 = (lt+273.16)**2
            a2 = (0.4156/(lt+273.16))
            a3 = 0.0020667

            evpf = (SP.exp((1137./a1)-a2-a3)-1)*1000; # equilibrium vapour pressure fractionation (per mil)

            do18e = (((1.+(evpf/1000.))*(1.+(tdf/1000.)+(((do18wvs/1000.)-(tdf/1000.))*(avpa/vpl))))-1.)*1000. # Craig-Gordon evaporation-site water (D\18Oe) (per mil wrt source water )
            do18l = (do18e*(1.-(SP.exp(-pe))))/pe # Bulk leaf water (D18OL) (per mil wrt source water)
            do18s = do18l+efc # Sucrose D18O (per mil wrt source water)
            do18c = do18l*(1.-(fwxm*feo))+efc # Cellulose D18O (per mil wrt source water)

            #convert back into ratios wrt SMOW:
            M = 1000.
            o18e = (((do18e/M+1.)*rsw)/rsmow-1.)*M
            o18l = (((do18l/M+1.)*rsw)/rsmow-1.)*M
            o18s = (((do18s/M+1.)*rsw)/rsmow-1.)*M
            o18c = (((do18c/M+1.)*rsw)/rsmow-1.)*M

            # correct to observed o18 given observed mean (c7)
            # and fraction (c8) water from meteoric sources
            o18cc = o18c
            o18cc = c8 * (o18c-(np.mean(o18c) - c7)) + (1 - c8) * c7
            self.dcell = o18cc

            # we also return dcell as a convenience if desired
            return self.dcell
