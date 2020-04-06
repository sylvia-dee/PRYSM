# Sylvia Dee
# PSM d18O Speleothem Calcite
# SENSOR SUB-MODEL
# Created 07/02/2014 <sdee@usc.edu>
# Modified 10/28/2014 <julieneg@usc.edu>

# import necessary modules
import numpy as np
import scipy.integrate as si
import scipy.special as sp
import psm.aux_functions.butter_lowpass_filter as bwf
# define utility functions for lowpass filtering

def adv_disp_transit(tau,tau0,Pe):
    #  advection-dispersion model [3]
    z0 = -0.5*np.sqrt(Pe*tau/tau0)   #(Eq 12a)
    zL = np.sqrt(Pe/tau*tau0) + z0   #(Eq 12b)
    h  = np.sqrt(1/(Pe*4*np.pi*tau/tau0))*np.exp(-0.25*Pe*tau/tau0)*(1 - np.exp(z0**2-zL**2)) + 0.25*(sp.erf(zL) - sp.erf(z0))  # Eq(11)
    return h

#  this is the sensor model
def speleo_sensor(t,d18O,T,dt,model='Adv-Disp',tau0=1.0,Pe=1.0):
    '''
    Speleothem Calcite [Sensor] Model
    Converts environmental signals to d18O calcite/dripwater using
    various models of karst aquifer recharge and calcification.

    INPUTS:
        dt        time step [years]                         (int or float, 1=1 year, 1./12. for months. etc.)
        t        time axis [years]                          (numpy array, length n)
        T        Average Annual Temperature     [K]         (numpy array, length n)
        d18O    delta-18-O (precip or soil water) [permil]  (numpy array, length n)
        NB: make sure these are precipitation weighted

    MODEL PARAMETERS
    ================
        model:  aquifer recharge model. possible values 'Well-Mixed'[1,2] or 'Adv-Disp' [3]
        tau0: mean transit time, [years]    (default = 0.5)
        Pe: Peclet number [non-dimensional] (default = 1.0) ('Adv-Disp' only)

     OUTPUTS:
         d18Oc  cave dripwater delta-18-O
         d18OK  delta-18-O of water after passage through karst
         h      transit time distribution

    REFERNCES:
        [1] Gelhar and Wilson: Ground Water Quality Modeling, Ground Water, 12(6):399--408
        [2] Partin et al., Geology, 2013, doi:10.1130/G34718.1  (SI)
        [3] Kirchner, J. W., X. Feng, and C. Neal (2001), Catchment-scale advection and dispersion as a mechanism for fractal scaling in stream
        [4] Wackerbarth et al. Modelling the d18O value of cave drip water and speleothem calcite, EPSL, 2010    doi:10.1016/j.epsl.2010.09.019
    '''
    #======================================================================
    # 1. Karst Model
    #======================================================================
    # ensure that tau0 is scaled correctly
    # (n.b. if in units of months, scale to years, *1/12)

    tau0=tau0*dt

    # Check to make sure mean residence time does not exceed 0.5 x [time series length]
    # It makes little sense to compute mixing with a long residence time with an extremely short timeseries.

    timeseries_len=len(d18O)#*dt
    if (tau0>=(0.5*timeseries_len)):
        print("Warning: mean residence time (tau0) is too close to half of the length of the timeseries.")



    # establish kernel length for transit time distribution
    # tau=np.arange(20.*tau0) {SDEE: revised 11/17/15}

    h_s=1E-4                        # ensure kernel approaches 0
    tau_s=-tau0*np.log(tau0*h_s)
    tau=np.arange(tau_s)

    # np.convolve method will freakout if the length of the kernel approaches that of the series.
    #return an error message and not perform the calculation: the effect of a kernel
    #of length L is roughly felt over 2*L +1, so 7 years in this case, which is the the length of the series.
    #Physically, it makes little sense to try and compute this in the first place,
    #so it would be better to explain to the user why this shouldn’t be attempted rather than return a meaningless estimate.

    if (len(tau)>=(0.5*timeseries_len)):
        print("Error: kernel length is too close to half of the length of the timeseries. This will result in spurious numerical effects: the effect of a kernel of length L is roughly felt over 2*L+1. Consider revising.")


    # define the transit time distribution h(tau) owing to several models
    if model == 'Well-Mixed':  # well-mixed reservoir model as in [1,2]
        h = (1./tau0) * np.exp(-tau/tau0)
    elif model == 'Adv-Disp':  #  advection-dispersion model [3]
        h = adv_disp_transit(tau,tau0,Pe)
        h[0] = adv_disp_transit(1e-3,tau0,Pe) # fix the pathological case at tau = 0
    else:
        print("Error: model " + model + " not found")

    # obtain normalization constant using Simpson's rule
    hint = si.simps(h,tau)
    # normalize transit time distribution
    h = h/hint

    # Normalize h to the min value of the timeseries
    #mindiff=min(np.abs(d18O-d18Om))

    # Remove mean to avoid boundary effects. (then put it back)
    d18Om = np.mean(d18O)

    # Shrink the exponential such that it doesn't add spurious variance to the resulting dripwater signal.
    hmax=max(np.abs(d18O-d18Om))
    h=h/hmax

    #convolve d18O input with the transit time distribution (Green's function of the problem)
    # this yields the d18O of aquifer water, after passage through the karst

    #d18OK = np.convolve((d18O-d18Om),h,mode='same') + d18Om
    d18OK = np.convolve((d18O-d18Om),h,mode='full') + d18Om

    # return to original shape
    #d18OK = d18OK[0:Ldata]
    #======================================================================
    #   2. Dripwater Model
    #======================================================================
    # NB: right now this is only temperature-Dependent Fractionation, but one could
    #    add any number of relevant processes here, as in [4]
    # NB: we assume yearly data, so T is the yearly average temperature for each time step.
    #   the "average temperature" used for the thermal fractionation calculation is
    #   somewhat arbitrarily the decadal average.
    #======================================================================
    # SDEE: use next part for yearly data/longer term speleos only.
    #avg_period = 10.0
    #Tm = np.mean(T)  # extract mean
    # lowpass filter the  temperature
    #Tl = bwf.filter(T-Tm,1.0/avg_period) + Tm
    #Tl = filter(T-Tm,1.0/avg_period) + Tm
    #======================================================================
    # apply thermal fractionation (Eq 11 in [4])
    Tl = T
    d18Oc = (d18OK + 1000)/1.03086*np.exp(2780/Tl**2 - 0.00289)-1000.0
#======================================================================
    return (d18Oc, d18OK, h, tau)
