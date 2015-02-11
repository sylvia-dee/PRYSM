# Sylvia Dee <sdee@usc.edu>
# PSM d18O Ice Core d18O
# Post-Processing (Spectral Analysis)
# Modified [sdee] 02/09/2015
# Script to Reproduce PRYSM Paper Figure 7
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


# Ice Cores (Quelccaya, Peru is the example in the paper figure).

d18Oi=np.load('Ice_Archive_d18O.npy')
IT=np.load('Ice_Temp.npy')
IP=np.load('Ice_Precip_Accum.npy')
Id18Op=np.load('Ice_Precip_d18O.npy')
#======================================================================================

# Compute Spectra for all variables

t=np.arange(1000,2005,1)
dt=1.0

# Run BAM SIMUL
#X=d18Oi.reshape(1005,1)
#tp, Xp, tmc =bam_simul_perturb(X,t,param=[0.01,0.01])
#np.save('QUEL_BAM_realiz.npy',Xp)

# Load Realizations

Xp=np.load('Ice_Age_Realizations.npy')

# Remove Mean for all variables

IT=IT-np.mean(IT)
IP=IP-np.mean(IP)
Id18Op=Id18Op-np.mean(Id18Op)
d18Oi=d18Oi-np.mean(d18Oi)
Xpm=Xp-np.mean(Xp, axis=0)

#
# Compute Spectra For all variables

ITf, ITpsd_mt, ITnu = tsa.multi_taper_psd(IT, Fs=1.0,adaptive=False, jackknife=False)
IPf, IPpsd_mt, IPnu = tsa.multi_taper_psd(IP, Fs=1.0,adaptive=False, jackknife=False)
Id18Opf, Id18Oppsd_mt, Id18Opnu = tsa.multi_taper_psd(Id18Op, Fs=1.0,adaptive=False, jackknife=False)
d18Oif, d18Oipsd_mt, d18Oinu = tsa.multi_taper_psd(d18Oi, Fs=1.0,adaptive=False, jackknife=False)

Xpf=np.zeros((len(ITf),len(Xpm[1])))
Xpsd_mt=np.zeros((len(ITf),len(Xpm[1])))
Xpnu=np.zeros((len(ITf),len(Xpm[1])))

for i in range(len(Xpm[1])):
    Xpf[:,i],Xpsd_mt[:,i],Xpnu[:,i]=tsa.multi_taper_psd(Xpm[:,i], Fs=1.0,adaptive=False, jackknife=False)


# Compute quantiles for spectra

from scipy.stats.mstats import mquantiles
q1=mquantiles(Xpsd_mt,prob=[0.025,0.975],axis=1)

# x axis for quantile of spectra:
q2=Xpf[:,0]
#======================================================================================

# Plot ice
#======================================================================================


# ENV

ax1=plt.subplot(4,2,2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.loglog(ITf,ITpsd_mt, label='T',color='red')
plt.loglog(IPf,IPpsd_mt, label='P',color='blue')
plt.title(r'ENVIRONMENT',fontsize=11, color='gray')
plt.xlabel(r'Period (Years)', fontsize=11)
plt.ylabel(r'PSD',fontsize=11)
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False)  
# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])
xtick=(1.0/pertick)
plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.minorticks_off()
plt.xlim([1./550.,0.4])
plt.ylim([1e-2, 1e1])
plt.legend(loc=3,fontsize=10,frameon=False)
#======================================================================================

ax2=plt.subplot(4,2,4)

plt.loglog(Id18Opf,Id18Oppsd_mt, label=r'Sensor $\delta^{18}O_{ICE}$',color='Indigo')

plt.title(r'SENSOR',fontsize=11, color='gray')
plt.xlabel(r'Period (Years)', fontsize=11)
plt.ylabel(r'PSD', fontsize=11)
ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False)  
# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.minorticks_off()
plt.xlim([1./550.,0.4])
plt.ylim([1e-2, 1e1])
plt.legend(loc=3,fontsize=11,frameon=False)
#======================================================================================
# ARCHIVE
ax3=plt.subplot(4,2,6)
plt.loglog(d18Oif,d18Oipsd_mt,label=r'$\delta^{18}O_{ICE}$ + Diffusion, Compaction',color='DodgerBlue',linewidth=1)

plt.title(r'ARCHIVE',fontsize=11, color='gray')
plt.xlabel(r'Period (Years)', fontsize=11)
plt.ylabel(r'PSD', fontsize=11)

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])
ax3.spines["top"].set_visible(False)  
ax3.spines["right"].set_visible(False)  
xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.minorticks_off()
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-4, 1])
plt.legend(loc=3,fontsize=11,frameon=False)

#======================================================================================
# OBS
ax4=plt.subplot(4,2,8)
#for i in range(len(Xpm[1])):
#    plt.loglog(Xpf[:,i],Xpsd_mt[:,i],color='gray')

plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations, CI',facecolor='gray',alpha=0.5)
plt.loglog(d18Oif,d18Oipsd_mt,label=r'$\delta^{18}O_{ICE}$ + Diffusion, Compaction',color='DodgerBlue',linewidth=1)


plt.title(r'OBSERVATION',fontsize=11, color='gray')
plt.xlabel(r'Period(Years)', fontsize=11)
plt.ylabel(r'PSD', fontsize=11)
ax4.spines["top"].set_visible(False)  
ax4.spines["right"].set_visible(False)  
# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.minorticks_off()
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-4, 1])
p1 = Rectangle((0, 0), 1, 1, fc="gray")
p2 = Rectangle((0, 0), 0.5, 0.5, fc="DodgerBlue")
legend([p1,p2],['1000 Age-Perturbed Realizations [2$\%$]',r'$\delta^{18}O_{ICE}$ + Diffusion, Compaction'],loc=3,fontsize=11,frameon=False)
#======================================================================================


#plt.suptitle(r'\textbf{MTM-Estimated Spectra: Ice Core PSM, Quelccaya}', fontsize=13)
plt.tight_layout()
plt.show()