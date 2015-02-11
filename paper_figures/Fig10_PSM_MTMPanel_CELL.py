# Sylvia Dee <sdee@usc.edu>
# PSM d18O Speleothem Calcite
# Post-Processing (Spectral Analysis)
# Modified [sdee] 02/09/2015
# Script to Reproduce PRYSM Paper Figure 10
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


cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Cellulose_Development/FigureVars/

# Trees (La Selva, Costa Rica)

TT=np.load('LaSelva_T.npy')
TP=np.load('LaSelva_P.npy')
Td18Op=np.load('LaSelva_d18OP.npy')
Td18Os=np.load('LaSelva_d18OS.npy')
TV=np.load('LaSelva_d18OV.npy')
TRH=np.load('LaSelva_RH.npy')
TC=np.load('LaSelva_d18Ocell.npy')
#======================================================================================

# Run Bam_SIMUL
#t=np.arange(1000,2005,1)
#X=TC.reshape(1005,1)
#tp, Xp, tmc =bam_simul_perturb(X,t,param=[0.02,0.02])
#np.save('LaSelva_AgePerturbed_Network.npy',Xp)
Xp=np.load('LaSelva_AgePerturbed_Network.npy')

# Compute Spectra for all variables

t=np.arange(1000,2005,1)
dt=1.0

# Remove Mean for all variables

TT=TT-np.mean(TT)
TP=TP-np.mean(TT)
Td18Op=Td18Op-np.mean(Td18Op)
TV=TV-np.mean(TV)
TRH=TRH-np.mean(TRH)
TC=TC-np.mean(TC)
TS=Td18Os-np.mean(Td18Os)

Xpm=Xp-np.mean(Xp, axis=0)

#
# Compute Spectra For all variables

TTf, TTpsd_mt, TTnu = tsa.multi_taper_psd(TT, Fs=1.0,adaptive=False, jackknife=False)
TPf, TPpsd_mt, TPnu = tsa.multi_taper_psd(TP, Fs=1.0,adaptive=False, jackknife=False)
Td18Opf, Td18Oppsd_mt, Td18Opnu = tsa.multi_taper_psd(Td18Op, Fs=1.0,adaptive=False, jackknife=False)
TVf, TVpsd_mt, TVnu = tsa.multi_taper_psd(TV, Fs=1.0,adaptive=False, jackknife=False)
TRHf, TRHpsd_mt, TRHnu = tsa.multi_taper_psd(TRH, Fs=1.0,adaptive=False, jackknife=False)
TCf, TCpsd_mt, TCnu = tsa.multi_taper_psd(TC, Fs=1.0,adaptive=False, jackknife=False)
TSf, TSpsd_mt, TSnu = tsa.multi_taper_psd(TS, Fs=1.0,adaptive=False, jackknife=False)

Xpf=np.zeros((len(TTf),len(Xpm[1])))
Xpsd_mt=np.zeros((len(TTf),len(Xpm[1])))
Xpnu=np.zeros((len(TTf),len(Xpm[1])))

for i in range(len(Xpm[1])):
    Xpf[:,i],Xpsd_mt[:,i],Xpnu[:,i]=tsa.multi_taper_psd(Xpm[:,i], Fs=1.0,adaptive=False, jackknife=False)


# Compute quantiles for spectra

from scipy.stats.mstats import mquantiles
q1=mquantiles(Xpsd_mt,prob=[0.025,0.975],axis=1)

# x axis for quantile of spectra:
q2=Xpf[:,0]
#======================================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot Trees
ax=plt.subplot(422)
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False) 
ax.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  


plt.loglog(TTf,TTpsd_mt, label='T',color='red')
plt.loglog(TPf,TPpsd_mt, label='P',color='blue')
plt.loglog(TRHf,TRHpsd_mt,label='RH',color='MediumVioletRed')

plt.title(r'ENVIRONMENT',fontsize=11, color='gray')
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
plt.ylim([1e-2, 1e2])
plt.legend(loc=3,ncol=3,mode='expand',fontsize=10,frameon=False)
#plt.legend(bbox_to_anchor=(0., -0.1, 1., -0.1), loc=3, ncol=4, mode="expand",
#======================================================================================

ax1=plt.subplot(424)
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False) 
plt.loglog(Td18Opf,Td18Oppsd_mt, label=r'$\delta^{18}O_{P}$',color='DarkCyan')
plt.loglog(TVf,TVpsd_mt, label=r'$\delta^{18}O_{V}$',color='MediumVioletRed')
plt.loglog(TSf,TSpsd_mt, label=r'$\delta^{18}O_{S}$',color='Sienna')

plt.title(r'WATER ISOTOPES',fontsize=11, color='gray')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
ax1.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax1.minorticks_off()
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-4, 1e1])
plt.legend(loc=3,ncol=3,mode='expand',fontsize=10,frameon=False)
#======================================================================================

ax2=plt.subplot(426)
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False) 
plt.loglog(TCf,TCpsd_mt,label=r'$\delta^{18}O_{C}$',color='DarkGreen')
plt.title(r'SENSOR',fontsize=11, color='gray')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
ax2.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax2.minorticks_off()
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-2, 1e0])
plt.legend(loc=3,fontsize=10,frameon=False)
plt.xticks(xtick, pertick_labels, size='small')

#======================================================================================


ax3=plt.subplot(428)
ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
# Observation purturbed age ensemble here.
plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations, CI',facecolor='gray',alpha=0.5)
plt.loglog(TCf,TCpsd_mt,label=r'$\delta^{18}O_{C}$',color='DarkGreen')
plt.xlabel(r'Frequency (Years)')
plt.ylabel(r'PSD')
ax3.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
ax3.minorticks_off()
plt.title(r'OBSERVATION',fontsize=11, color='gray')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.xticks(xtick, pertick_labels, size='small')

plt.ylim([1e-2, 1e0])
p2 = Rectangle((0, 0), 1, 1, fc="gray")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="DarkGreen")
legend([p1,p2],[r'$\delta^{18}O_{C}$','1000 Age-Perturbed Realizations [4$\%$], CI'],loc=3,fontsize=10,frameon=False)

#======================================================================================


#plt.suptitle(r'Tree Ring $\delta^{18}O_{C}$, La Selva',fontsize=14, color='gray')
plt.tight_layout()
plt.show()
