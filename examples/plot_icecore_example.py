# ======================================================================
# PLOTTING EXAMPLES
# ======================================================================
import numpy as np
import matplotlib.pyplot as plt

###Load data
outdir='./results/'
ice_Xn=np.load(outdir+"ice_Xn.npy")
ice_diffused=np.load(outdir+"ice_diffused.npy")
z=np.load(outdir+"ice_depth.npy")
time_d=np.load(outdir+"ice_time_d.npy")
###

# Example 1: Diffusion Profiles

plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.rcParams['font.family']='serif'

diffs =((np.diff(z)/np.diff(time_d)))
diffs=diffs[0:-1]
# row and column sharing
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.plot(z[0:-1],D*1e4)
ax2.plot(z[0:-2],sig*100.0,color='c')
ax3.plot(z,time_d,color='r')
ax4.plot(z[0:-2],sig/diffs,color='g')
# Set common labels
ax1.set_ylabel(r'diffusivity ($cm^2/s$)')
ax1.set_xlabel('Depth(m)')
ax2.set_ylabel('Diffusion length (cm)')
ax2.set_xlabel('Depth(m)')
ax3.set_ylabel('Age (yrs)')
ax3.set_xlabel('Depth (m)')
ax4.set_ylabel('Diffusion Length/Annual Layer Thickness')
ax4.set_xlabel('Depth (m)')
plt.suptitle('Quelccaya Diffusion Lengths')
plt.show()

#======================================================================

# Example 2: Plot age realizations with depth (THIS IS THE SAVED XP)
# Change values to fit your data in plotting examples.

bam=Xp
bam.shape

depths=depth_horizons
depths=np.flipud(depths)

quel=ice_diffused

for i in range(len(bam[1])):
    plt.plot(depths,bam[:,i],color='0.75')

plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.rcParams['font.family']='serif'

plt.plot(depths,quel,color='Blue')
#plt.legend(loc=3,fontsize=14,frameon=False)
#plt.xlim([1000,2010])
plt.ylim([-12,-10.5])
plt.xlabel(r'Depth in Ice Core (m)')
plt.ylabel(r'Modeled $\delta^{18}O_{ICE}$')
#plt.gca().invert_yaxis()
plt.title(r'\texttt{BAM} Simulated Chronologies, 2$\%$ Symmetric Dating Errors, Quelccaya Ice Cap, Peru')
plt.show()
