# Sylvia Dee <sdee@usc.edu>
# PRYSM
# Post-Processing and Figures
# Figure 5: Quelccaya Data-Model Comparison
# Script to reproduce paper results
# Modified [sdee] 01/30/2015
#===================================================================

# NOTE TO USER: THIS SCRIPT IS DESIGNED TO BE RUN IN IPYTHON CONSOLE
# highlight, copy all text, >>%paste

# Plot Quelccaya real data and simulated data. How do they compare?

#===================================================================
from pydap.client import open_url  
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


from __future__ import unicode_literals

import nitime.algorithms as tsa
import nitime.utils as utils
from nitime.viz import winspect
from nitime.viz import plot_spectral_estimate
#===============================================================================
# Load Datasets

#cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Ice_Development/TesterVars

# set up time axis
years=np.arange(1000,2005,1)

# LOAD SPEEDY DATA

quel=np.load('SPEEDY_quel_d18O_Ice.npy')
p=np.load('SPEEDY_quel_d18O_precip.npy')

# LOAD ECHAM5-wiso DATA
ECHAM_p=np.load('ECHAM_quel_d18O_precip_orig.npy')
ECHAM_piso=np.load('ECHAM_quel_d18O_precip.npy')
ECHAM_ice =np.load('ECHAM_quel_d18O_Ice.npy')

# LOAD MEASURED DATA (Thompson et al. 2013)
year=np.load('DATA_time.npy')
quel2=np.load('DATA_d18O_ice.npy')

#===============================================================================

# Compute spectra for all three datasets: SPEEDY-IER simulated ice, ECHAM, Real

#ECHAM = ECHAM5 wiso simulated for 1871-2011
#quel2 = real data
#quel  = speedy ice
#p     = d18O p for SPEEDY

# Remove Mean for all variables

echam_p=ECHAM_p-np.mean(ECHAM_p)
echam_piso=ECHAM_piso-np.mean(ECHAM_piso)
echam_ice=ECHAM_ice-np.mean(ECHAM_ice)
quel_data=quel2-np.mean(quel2)
quel_speedy=quel-np.mean(quel)
p_speedy=p-np.mean(p)

#
# Compute Spectra For all variables

ef, epsd_mt, enu = tsa.multi_taper_psd(echam_p, Fs=1.0,adaptive=False, jackknife=False)
epf, eppsd_mt, epnu = tsa.multi_taper_psd(echam_piso, Fs=1.0,adaptive=False, jackknife=False)
eif, eipsd_mt, einu = tsa.multi_taper_psd(echam_ice, Fs=1.0,adaptive=False, jackknife=False)

queldf, queldpsd_mt, queldnu = tsa.multi_taper_psd(quel_data, Fs=1.0,adaptive=False, jackknife=False)
quelsf, quelspsd_mt, quelsnu = tsa.multi_taper_psd(quel_speedy, Fs=1.0,adaptive=False, jackknife=False)
pf, ppsd_mt, pnu = tsa.multi_taper_psd(p, Fs=1.0,adaptive=False, jackknife=False)


# plot together
#===============================================================================

time1=np.arange(1000,2005,1)
time2=np.arange(1871,2012,1)
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.rcParams['font.family']='serif'

ax1=plt.subplot(2,1,1)
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False) 
ax1.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.plot(year,quel2,label=r'Measured $\delta^{18}O_{ICE}$ Data (Thompson 2013), $\sigma=2.1$',color='Navy',linewidth=1)

plt.plot(time1,p,label=r'SPEEDY-IER + PSM-Simulated $\delta^{18}O_{P}$, $\sigma=0.07$',color='Indigo',linewidth=1)
plt.plot(time1,quel,label=r'SPEEDY-IER + PSM-Simulated $\delta^{18}O_{ICE}$',color='violet',linestyle='--',linewidth=1)
plt.plot(time2,ECHAM_p,label=r'ECHAM5-wiso-Simulated $\delta^{18}O_{P}$, $\sigma=0.9$',color='DodgerBlue',linewidth=1)
plt.plot(time2,ECHAM_piso,label=r'ECHAM5-wiso-Simulated $\delta^{18}O_{P}$ + Altitude Correction',color='DeepSkyBlue',linewidth=1)
plt.plot(time2,ECHAM_ice,label=r'ECHAM5-wiso-Simulated $\delta^{18}O_{ICE}$, $\sigma=0.9$',linestyle='--',color='blue',linewidth=1)


plt.title('Timeseries Data',fontsize=11, color='gray')
ax1.minorticks_off()
plt.legend(loc=3,ncol=1,fontsize=11,frameon=False)
plt.xlim([1000,2011])
plt.ylim([-40,-7])
plt.ylabel(r'$\delta^{18}O$ [(\u2030)]')
plt.xlabel('Year')

#===============================================================================

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot 
ax=plt.subplot(2,1,2)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False) 
ax.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  

plt.loglog(queldf,queldpsd_mt, label=r'Measured $\delta^{18}O_{ICE}$ Data (Thompson 2013)',color='Navy',linewidth=1)

plt.loglog(pf,ppsd_mt, label=r'SPEEDY $\delta^{18}O_P$',color='Indigo',linewidth=1)
plt.loglog(quelsf,quelspsd_mt,label=r'SPEEDY $\delta^{18}O_{ICE}$ (+ Diffusion)',color='Indigo',linestyle='--',linewidth=1)

plt.loglog(epf,eppsd_mt, label=r'ECHAM5-wiso $\delta^{18}O_P$',color='DeepSkyBlue',linewidth=1)
plt.loglog(eif,eipsd_mt, label=r'ECHAM5-wiso $\delta^{18}O_{ICE}$ (+ Diffusion)',color='DeepSkyBlue',linestyle='--',linewidth=1)


plt.title(r'MTM-Generated Spectra',fontsize=11, color='gray')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
ax.minorticks_off()
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-5, 1e2])
plt.legend(loc=3,ncol=1,fontsize=11,frameon=False)

plt.suptitle(r'Model-Data Comparison for $\delta^{18}O_{ICE}$, Quelccaya Ice Cap, Peru')
plt.show()