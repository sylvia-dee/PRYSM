# Sylvia Dee 
# PSM d18O Ice Core
# ARCHIVE SUB-MODEL
# Modified 01/26/2015 <sdee@usc.edu>

# =======================================================================================

def diffusivity(rho,T=250,P=0.7,rho_d=822,b=0.25):

    '''
    DOCSTRING: Function 'diffusivity'
    Description: Calculates diffusivity (in m^2/s) as a function of density.
    
    Inputs:
    P: Ambient Pressure in Atm
    T: Temperature in K
    rho: density profile (kg/m^3)
    rho_d: 822 kg/m^2 [default], density at which ice becomes impermeable to diffusion
    
    Defaults are available for all but rho, so only one argument need be entered. 
    
    Note values for diffusivity in air:
    
    D16 = 2.1e-5*(T/273.15)^1.94*1/P
    D18 = D16/1.0285
    D2 = D16/1.0251
    D17 = D16/((D16/D18)^0.518)
    
    Reference: Johnsen et al. (2000): Diffusion of Stable isotopes in polar firn and ice:
    the isotope effect in firn diffusion
        
    '''
    import numpy as np
    import scipy
    from scipy import integrate
    import matplotlib.pyplot as plt

    # Set Constants
    
    R = 8.314478                                                # Gas constant
    m = 18.02e-3                                                # molar weight of water (in kg)
    alpha18 = np.exp(11.839/T-28.224e-3)                    # ice-vapor fractionation for oxygen 18
    p=np.exp(9.5504+3.53*np.log(T)-5723.265/T-0.0073*T)     # saturation vapor pressure
    Po = 1.                                                 # reference pressure, atmospheres
    ppa=3.454e12*np.exp(-6133/T)
    rho_i = 920.# kg/m^3, density of solid ice
    
    # Set diffusivity in air (units of m^2/s)
    
    Da = 2.1e-5*np.power((T/273.15),1.94)*(Po/P)        
    Dai = Da/1.0285

    # Calculate Tortuosity 
    
    invtau=np.zeros(len(rho))
    for i in range(len(rho)):
        if rho[i]<=rho_i/np.sqrt(b):
            invtau[i]=1-b*np.power((rho[i]/rho_d),2)
        else:
            invtau[i]=0


    D=m*p*invtau*Dai*(1/rho-1/rho_d)/(R*T*alpha18)

    return D
    
# =======================================================================================
# =======================================================================================

def icecore_diffuse(d18O,b,time,T,P,depth,depth_horizons,dz,drho):

    '''
        DOCSTRING: Function 'icecore_diffuse'
        DESCRIPTION: accounts for diffusion and compaction in the firn.

        Inputs:
            d18O: ice core isotope ratio, output from sensor Model (permil)
            b   : average accumulation rate at site (m/year)
            time: calendar years of record
            T   : average temperature at site (K)
            P   : average sea level pressure at site (atm)
            depth: total depth of core
            depth_horizons: accumulation by year in ice core--depth horizons (moving downwards in core) in meters
            dz: step in depth (default = min(depth_horizons)/10.) 
            drho: step in density (default=0.5 kg/m^3)

        Diffusion is computed using a convolution (Gaussian smoothing).

        Functionality: Calculates diffusion length as a function of density
        in firn, and given vectors of time-depth and density-depth 
        Also expects Pressure in atm, T in K, rho and rho_ice in 
        kg/m^3, but defaults on these, so only the first three arguments
        must be entered.

        "Time" is really "age" (increasing down core) and is given in years.
        z is depth in meters
        rho should be in kg/m^3
        the vectors rho, time, and z should all correponding with one another
    '''
# ======================================================================
# A.0: Initialization
# ====================================================================== 
    import numpy as np
    import scipy
    from scipy import integrate
    from scipy import ndimage
    from scipy import stats
    import matplotlib.pyplot as plt
    from math import pi, sqrt, exp

    # Set constants (note units all in m)

    R=8.314478# Gas constant
    m = 18.02e-3# molar weight of water (in kg)
    rho_s = 300.#kg/m^3, surface density
    rho_d = 822.#kg/m^2, density at which ice becomes impermeable to diffusion
    rho_i = 920.# kg/m^3, density of solid ice

# ======================================================================
# A.1: Compaction Model
# ======================================================================

# Set Density Profile
    z = np.arange(0,depth,dz)+dz# linear depth scale
    rho = np.linspace(rho_s,rho_d,len(z))
    rho = rho[0:len(z)]
    time_d = np.cumsum(dz/b*rho/rho_i)# timescale accounting for densification
    ts = time_d*365.25*24*3600# convert time in years to ts in seconds

    # Integrate diffusivity along the density gradient to obtain diffusion length

    drho = np.diff(rho)
    dtdrho = np.diff(ts)/np.diff(rho)
    D=diffusivity(rho,T,P,rho_d,b)
    D=D[0:-1]
    rho = rho[0:-1]
    diffs =((np.diff(z)/np.diff(time_d)))
    diffs=diffs[0:-1]

    # Integration using the trapezoidal method

    sigma_sqrd_dummy = 2*np.power(rho,2)*dtdrho*D*drho
    sigma_sqrd = integrate.cumtrapz(sigma_sqrd_dummy)
    rho=rho[0:-1]
    sigma = np.sqrt(1/np.power(rho,2)*sigma_sqrd)

# ====================================================================== 
# A.2. Diffusion Profile 
# ====================================================================== 

    # Load water isotope series

    del18 = np.flipud(d18O)# NOTE YOU MIGHT NOT NEED FLIP UD here. Our data goes forward in time.
    depth_horizons = depth_horizons

    # interpolate over depths to get an array of dz values corresponding to isotope values
    # for convolution/diffusion

    years_rev=np.flipud(time)# again, may not need to invert time axis here.

    iso_interp=np.interp(z,depth_horizons,del18)
    time_interp=np.interp(z,depth_horizons,years_rev)

    diffused_final=np.zeros(len(iso_interp))

#======================================================================
    # Establish Kernel for convolution of timeseries with a Gaussian
#======================================================================
    # start out with a very wide kernel and then clip to ensure accurate sampling of Gaussian. 
    # return a warning if the kernel length is approaching 1/2 that of the timeseries. This will result in spurious numerical effects.
    
    zp=np.arange(-100,100,dz)
    if (len(zp) >= 0.5*len(z)):
        print "Warning: convolution kernal length (zp) is approaching that of half the length of timeseries. Kernal being clipped."
        bound=0.249*len(z)*dz
        zp=np.arange(-bound,bound,dz)

#======================================================================
    # Do the entire calculation on the isotope time series

    for i in range(len(sigma)):
        part1=(1./(sigma[i]*np.sqrt(2.*np.pi)))
        part2=scipy.exp(-zp**2/(2*sigma[i]**2))
        G=part1*part2
        #remove mean and then put back
        rm = np.mean(iso_interp)
        cdel = iso_interp-rm
        diffused = np.convolve(G,cdel,mode='same')*dz
        diffused = diffused + rm
        diffused_final[i]=diffused[i]

# take off the first few and last few points used in convolution
    diffused_timeseries = diffused_final[0:-3]

# Now we need to pack our data back into single year data units based on the depths and year interpolated data
    final_iso=np.interp(depth_horizons,z[0:-3],diffused_timeseries)
    ice_diffused=final_iso
# ====================================================================== 
    return z, sigma, D, time_d, diffs, ice_diffused


