
# Sylvia Dee <sylvia_dee@brown.edu>
# PRYSM
# PSM for Lacustrine Sedimentary Archives
# ARCHIVE MODEL: Compaction of Sediments based on porosity
# Function 'porosity'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
# Modified 04/2/2018 <sylvia@ig.utexas.edu>

#======================================================================================

def porosity(depth_profile):

    
    import numpy as np    
    # assuming no excess pore pressure (Pe = 0) and that grain density is constant for quartz:
    z=np.linspace(0,depth_profile,num=100) #depth, meters
    ps=2650.     # kg/m^3 (grain density of quartz, valid for most sands/seds)
    pw = 1000.   # kg/m^3
    g = 9.8     # m/s^2 (gravity)
    phi=np.zeros(len(z))
    phi_0=0.95  # surface porosity, set by user (for reference, shale ~ 0.6, 0.3 to 0.5 for sands, 0.2-0.5 for silts and shales)

    k1 = (1-phi_0)/phi_0
    c = 3.68E-8 
    for i in range(len(phi)):
        phi[i] = (np.exp(-c*g*(ps-pw)*z[i]))/(np.exp(-c*g*(ps-pw)*z[i])+k1)

    return phi, z
