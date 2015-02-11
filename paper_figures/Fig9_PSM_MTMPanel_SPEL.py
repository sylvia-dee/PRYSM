# Sylvia Dee <sdee@usc.edu>
# PSM d18O Speleothem Calcite
# Post-Processing (Spectral Analysis)
# Modified [sdee] 02/09/2015
# Script to Reproduce PRYSM Paper Figure 9
#==============================================================
# E0. Initialization
#==============================================================

from pydap.client import open_url  
from pylab import *

from math import exp

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.stats.distributions as dist

import nitime.algorithms as tsa
import nitime.utils as utils
from nitime.viz import winspect
from nitime.viz import plot_spectral_estimate
#
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Load Data for Speleothems (Hidden Cave) 
cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Speleothem_Development/Figure_Vars/
ST=np.load('Cave_Temp.npy')
SP=np.load('Cave_Precip.npy')
Sd18Op=np.load('Cave_d18OP.npy')
Sd18Os=np.load('Cave_d18OS.npy')
Sd18Oc=np.load('Cave_d18O_Calcite.npy')
#======================================================================================

# Compute Spectra for all variables

t=np.arange(1000,2005,1)
dt=1.0

# Load Age-perturbed data:
Bchron=np.load('Cave_Age_Realizations.npy')
# Remove Mean for all variables
Xpm=Bchron-np.mean(Bchron, axis=0)
#
ST=ST-np.mean(ST)
SP=SP-np.mean(SP)
Sd18Op=Sd18Op-np.mean(Sd18Op)
Sd18Os=Sd18Os-np.mean(Sd18Os)
Sd18Oc=Sd18Oc-np.mean(Sd18Oc)
# Compute Spectra For all variables

STf, STpsd_mt, STnu = tsa.multi_taper_psd(ST, Fs=1.0,adaptive=False, jackknife=False)
SPf, SPpsd_mt, SPnu = tsa.multi_taper_psd(SP, Fs=1.0,adaptive=False, jackknife=False)
Sd18Opf, Sd18Oppsd_mt, Sd18Opnu = tsa.multi_taper_psd(Sd18Op, Fs=1.0,adaptive=False, jackknife=False)
Sd18Osf, Sd18Ospsd_mt, Sd18Osnu = tsa.multi_taper_psd(Sd18Os, Fs=1.0,adaptive=False, jackknife=False)
Sd18Ocf, Sd18Ocpsd_mt, Sd18Osnu = tsa.multi_taper_psd(Sd18Oc, Fs=1.0,adaptive=False, jackknife=False)

Xpf=np.zeros((len(STf),len(Xpm[1])))
Xpsd_mt=np.zeros((len(STf),len(Xpm[1])))
Xpnu=np.zeros((len(STf),len(Xpm[1])))

for i in range(len(Xpm[1])):
    Xpf[:,i],Xpsd_mt[:,i],Xpnu[:,i]=tsa.multi_taper_psd(Xpm[:,i], Fs=1.0,adaptive=False, jackknife=False)


# Compute quantiles for spectra

from scipy.stats.mstats import mquantiles
q1=mquantiles(Xpsd_mt,prob=[0.025,0.975],axis=1)

# x axis for quantile of spectra:
q2=Xpf[:,0]

#======================================================================================

# Plot Speleo
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax=plt.subplot(422)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False) 
ax.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.loglog(STf,STpsd_mt, label='T',color='red')
plt.loglog(SPf,SPpsd_mt, label='P',color='blue')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
plt.title(r'ENVIRONMENT',fontsize=11, color='gray')
ax.minorticks_off()

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-3, 1e1])
plt.legend(loc=3,ncol=2,fontsize=11,frameon=False)
#======================================================================================

ax2=plt.subplot(424)
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False) 
plt.loglog(Sd18Opf,Sd18Oppsd_mt, label=r'Sensor $\delta^{18}O_{P}$',color='MidnightBlue')
plt.title(r'SENSOR',fontsize=11, color='gray')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-2, 1e0])
plt.legend(loc=3,fontsize=11,frameon=False)
ax2.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax2.minorticks_off()
#======================================================================================

ax3=plt.subplot(426)
ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False) 
plt.loglog(Sd18Ocf,Sd18Ocpsd_mt,label=r'$\delta^{18}O_{C}$',color='Maroon')
plt.xticks(xtick, pertick_labels, size='small')
plt.xlim([1./550.,0.4])
plt.title(r'ARCHIVE',fontsize=11, color='gray')
plt.grid('on',axis='y',color='DimGray')

plt.ylim([1e-1, 1e1])
ax3.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax3.minorticks_off()

plt.legend(loc=3,fontsize=11,frameon=False)
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')

#======================================================================================

ax4=plt.subplot(428)
plt.loglog(Sd18Ocf,Sd18Ocpsd_mt,label=r'$\delta^{18}O_{C}$',color='Maroon')
plt.fill_between(q2,q1[:,0],q1[:,1],label=r'1000 \texttt{Bchron} Realizations, CI',facecolor='gray',alpha=0.5)

ax4.spines["top"].set_visible(False)  
ax4.spines["right"].set_visible(False) 
# BCRHON
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
ax4.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax4.minorticks_off()
plt.title(r'OBSERVATION',fontsize=11, color='gray')
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-1, 1e1])
p2 = Rectangle((0, 0), 1, 1, fc="gray")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="Maroon")
legend([p1,p2],[r'$\delta^{18}O_{C}$',r'1000 \texttt{Bchron} Realizations'],loc=3,fontsize=11,frameon=False)

#plt.suptitle(r'Speleothem $\delta^{18}O_{C}$, Borneo',fontsize=14, color='gray')
plt.tight_layout()
plt.show()
