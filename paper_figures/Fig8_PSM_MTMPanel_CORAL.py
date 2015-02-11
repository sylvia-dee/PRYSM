# Sylvia Dee <sdee@usc.edu>
# PSM d18O Coral Aragonite
# Post-Processing (Spectral Analysis)
# Plot Spectra for a coral site showing sub-model tranforms
# Modified [sdee] 02/09/2015
# Script to Reproduce PRYSM Paper Figure 8
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
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


# Corals (Palmyra)

cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Coral_Development/Figure_Vars/
SST=np.load('Palmyra_SST.npy')
SSS=np.load('Palmyra_SSS.npy')
Cd18O=np.load('Palmyra_d18O.npy')
#======================================================================================

# Compute Spectra for all variables

t=np.arange(1000,2005,1)
dt=1.0

# Bam Simul
#X=Cd18O.reshape(1005,1)
#tp, Xp, tmc =bam_simul_perturb(X,t,param=[0.025,0.025])

#======================================================================================

Xp=np.load('Palmyra_Age_Perturbed.npy')
# Remove Mean for all variables

CST=SST-np.mean(SST)
CSS=SSS-np.mean(SSS)
Cd18O=Cd18O-np.mean(Cd18O)
Xpm=Xp-np.mean(Xp, axis=0)

# Compute Spectra For all variables

CSTf, CSTpsd_mt, CSTnu = tsa.multi_taper_psd(CST, Fs=1.0,adaptive=False, jackknife=False)
CSSf, CSSpsd_mt, CSSnu = tsa.multi_taper_psd(CSS, Fs=1.0,adaptive=False, jackknife=False)
Cd18Of, Cd18Opsd_mt, Cd18Onu = tsa.multi_taper_psd(Cd18O, Fs=1.0,adaptive=False, jackknife=False)


Xpf=np.zeros((len(CSTf),len(Xpm[1])))
Xpsd_mt=np.zeros((len(CSTf),len(Xpm[1])))
Xpnu=np.zeros((len(CSTf),len(Xpm[1])))

for i in range(len(Xpm[1])):
    Xpf[:,i],Xpsd_mt[:,i],Xpnu[:,i]=tsa.multi_taper_psd(Xpm[:,i], Fs=1.0,adaptive=False, jackknife=False)


# Compute quantiles for spectra

from scipy.stats.mstats import mquantiles
q1=mquantiles(Xpsd_mt,prob=[0.025,0.975],axis=1)

# x axis for quantile of spectra:
q2=Xpf[:,0]
#======================================================================================
# Plot Corals
#======================================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax=plt.subplot(322)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)  
ax.loglog(CSTf,CSTpsd_mt, label='NINO3.4 SST',color='red')
ax.loglog(CSSf,CSSpsd_mt, label='SSS (PSU)',color='LightSeaGreen')
ax.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax.minorticks_off()
ax.set_title(r'ENVIRONMENT',fontsize=14, color='gray')
ax.set_xlabel(r'Period (Years)')
ax.set_ylabel(r'PSD')

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-3, 1e1])
plt.legend(loc=3,fontsize=11,frameon=False)

#======================================================================================

ax2=plt.subplot(324)
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False)  
plt.loglog(Cd18Of,Cd18Opsd_mt, label=r'Coral $\delta^{18}O_{C}$',color='DarkOrange')
plt.title(r'SENSOR',fontsize=14, color='gray')
plt.xlabel(r'Period (Years)')
plt.ylabel(r'PSD')
ax2.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax2.minorticks_off()
# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-3, 1])
ax2.legend(loc=3,fontsize=11,frameon=False)
#======================================================================================


ax3=plt.subplot(326)
ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
# Observation purturbed age ensemble here.
plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations (4%), CI',facecolor='gray',alpha=0.5)
plt.loglog(Cd18Of,Cd18Opsd_mt, label=r'Coral $\delta^{18}O_{C}$',color='DarkOrange')
ax3.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax3.minorticks_off()
plt.title(r'OBSERVATION',fontsize=14, color='gray')
plt.xlabel(r'Period (Years)')
plt.ylabel(r'PSD')

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-3, 1])
p2 = Rectangle((0, 0), 1, 1, fc="gray")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="DarkOrange")
legend([p1,p2],[r'Coral $\delta^{18}O_{C}$','1000 Age-Perturbed Realizations [4$\%$], CI'],loc=3,fontsize=11,frameon=False)

#======================================================================================

#plt.suptitle(r'\textbf{Sub-Model MTM-Estimated Spectra for Coral PSM, Palmyra Island}', fontsize=16)
plt.tight_layout()
plt.show()