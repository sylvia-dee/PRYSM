# Sylvia Dee 
# PSM d18O Ice Core
# ARCHIVE SUB-MODEL
# Modified 01/26/2015 <sdee@usc.edu>
# Modified 11/18/2015 <sylvia_dee@brown.edu>

# =======================================================================================

def diffusivity(rho,T=250,P=0.9,rho_d=822,b=1.3):
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
            #invtau[i]=1.-1.3*np.power((rho[i]/rho_d),2)
            invtau[i]=1.-1.3*np.power((rho[i]/rho_i),2)
        else:
            invtau[i]=0.


    D=m*p*invtau*Dai*(1/rho-1/rho_d)/(R*T*alpha18)

    return D

# =======================================================================================
# =======================================================================================

def densification(Tavg,bdot,rhos,z):#,model='hljohnsen'):

    '''## Calculates steady state snow/firn depth density profiles using Herron-Langway type models.
        
        # [rho,zieq]=densification(Tavg,bdot,rhos,z,model)
        
        # Herron-Langway type models. (Arthern et al. 2010 formulation).
        
        # INPUT:
        # Tavg: 10m temperature in celcius ## CELCIUS!
        # bdot: accumulation rate in mwe/yr or (kg/m2/yr)
        # rhos: surface density in kg/m3
        # z: depth in true_metres 
        
        # model can be: {'HLJohnsen' 'HerronLangway' 'LiZwally' 'Helsen' 'NabarroHerring'}
        # default is herronlangway. (The other models are tuned for non-stationary modelling (Read Arthern et al.2010 before applying in steady state).
        
        # OUTPUT: 
        # rho: density (kg/m3) for all z-values.
        # zieq: ice equivalent depth for all z-values.
        # t: age for all z-values (only taking densification into account.)
        
        # Example usage:
        # z=0:300
        # [rho,zieq,t]=densitymodel(-31.5,177,340,z,'HerronLangway')
        # plot(z,rho)
        # Aslak Grinsted, University of Copenhagen 2010
        # Adapted Sylvia Dee, Brown University, 2017
    '''
    import numpy as np
    import scipy
    from scipy import integrate
    import matplotlib.pyplot as plt
    
    rhoi=920. 
    rhoc=550.
    #rhoc=822.
    rhow=1000. 
    rhos = 340.
    R=8.314
    #tmodel=model

    #Tavg=248.
    #bdot=0.1
    #Herron-Langway with Johnsen et al 2000 corrections.
    #Small corrections to HL model which are not in Arthern et al. 2010 

    c0=0.85*11*(bdot/rhow)*np.exp(-10160./(R*Tavg))
    c1=1.15*575*np.sqrt(bdot/rhow)*np.exp(-21400./(R*Tavg))

    k0=c0/bdot ##~g4
    k1=c1/bdot

    #critical depth at which rho=rhoc
    zc=(np.log(rhoc/(rhoi-rhoc))-np.log(rhos/(rhoi-rhos)))/(k0*rhoi) #g6

    ix=z<=zc #find the z's above and below zc
    upix=np.where(ix) #indices above zc
    dnix=np.where(~ix) #indices below zc

    q=np.zeros((z.shape)) #pre-allocate some space for q, rho
    rho=np.zeros((z.shape))

    #test to ensure that this will not blow up numerically if you have a very very long core.
    # manually set all super deep layers to solid ice (rhoi=920)
    NUM=k1*rhoi*(z-zc)+np.log(rhoc/(rhoi-rhoc))

    numerical= np.where(NUM<=100.0)
    blowup=np.where(NUM>100.0)

    q[dnix]=np.exp(k1*rhoi*(z[dnix]-zc)+np.log(rhoc/(rhoi-rhoc))) #g7
    q[upix]=np.exp(k0*rhoi*z[upix]+np.log(rhos/(rhoi-rhos))) #g7

    rho[numerical]=q*rhoi/(1+q) #[g8]
    rho[blowup]=rhoi

    #only calculate this if you want zieq
    tc=(np.log(rhoi-rhos)-np.log(rhoi-rhoc))/c0 #age at rho=rhoc [g17]
    t=np.zeros((z.shape)) #pre allocate a vector for age as a function of z
    t[upix]=(np.log(rhoi-rhos)-np.log(rhoi-rho[upix]))/c0 # [g16] above zc
    t[dnix]=(np.log(rhoi-rhoc)-np.log(rhoi+0.0001-rho[dnix]))/c1+tc # [g16] below zc
    tdiff=np.diff(t)

    # make sure time keeps increasing even after we reach the critical depth.
    if np.any(tdiff==0.00):
        for i in range(len(t)):
            inflection=np.where(tdiff==0.0)
            lineardepth_change=t[inflection][0]

            for i in range(len(t)):
                if t[i]<lineardepth_change:
                    t[i]=t[i]
                elif t[i]>lineardepth_change:
                    t[i]=t[i-1]+1E-5

    zieq=t*bdot/rhoi #[g15]

    rho=np.array(rho)
    return rho,zieq,t

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
    import matplotlib.pyplot as plt

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
    #rho = np.linspace(rho_s,rho_d,len(z))

    # CALL DENSIFICATION
    rho,zieq,t=densification(T,b,rho_s,z)
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

    # IMPORTANT: once the ice reaches crtiical density (solid ice), there will no longer
    # be any diffusion. There is also numerical instability at that point. Set Sigma=1E-13 for all
    # points below that threshold.

    # Set to 915 to be safe.
    solidice=np.where(rho>=rho_d-5.0)
    diffusion=np.where(rho<rho_d-5.0)

    sigma_sqrd_dummy = 2*np.power(rho,2)*dtdrho*D*drho
    sigma_sqrd= integrate.cumtrapz(sigma_sqrd_dummy)
    rho=rho[0:-1]
    sigma=np.zeros((len(rho)+1))
    sigma[diffusion] = np.sqrt(1/np.power(rho,2)*sigma_sqrd)
    #sigma[solidice]=np.nanmax(sigma) #max diffusion length in base of core // set in a better way. max(sigma)
    sigma[solidice]=sigma[diffusion][-1]
    sigma=sigma[0:-1]
    #sigma=sigma*0.5
    # ====================================================================== 
    # A.2. Diffusion Profile 
    # ====================================================================== 

    # Load water isotope series

    del18 = np.flipud(d18O)# NOTE YOU MIGHT NOT NEED FLIP UD here. Our data goes forward in time.
    depth_horizons = depth_horizons

    # interpolate over depths to get an array of dz values corresponding to isotope values
    # for convolution/diffusion

    years_rev=np.flipud(time)# again, may not need to invert time axis here.
    z=np.reshape(z, len(z))
    del18=np.reshape(del18,len(del18))
    iso_interp=np.interp(z,depth_horizons,del18)
    time_interp=np.interp(z,depth_horizons,years_rev)

    diffused_final=np.zeros(len(iso_interp))

    #======================================================================

    # return a warning if the kernel length is approaching 1/2 that of the timeseries. This will result in spurious numerical effects.

    zp=np.arange(-100,100,dz)
    if (len(zp) >= 0.5*len(z)):
        print("Warning: convolution kernal length (zp) is approaching that of half the length of timeseries. Kernal being clipped.")
        bound=0.20*len(z)*dz
        zp=np.arange(-bound,bound,dz)

    #======================================================================

    sigma_dummy=np.tile(0.08,len(sigma))

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
    # ===================================================
    return z, sigma, D, time_d, diffs, ice_diffused, rho, zieq

