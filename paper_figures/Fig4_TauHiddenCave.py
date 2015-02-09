# Sylvia Dee <sdee@usc.edu>
# PRYSM
# Post-Processing and Figures
# Figure 4: Karst Processes
# Script to reproduce paper results
# Modified [sdee] 01/30/2015
#===================================================================

# NOTE TO USER: THIS SCRIPT IS DESIGNED TO BE RUN IN IPYTHON CONSOLE
# highlight, copy all text, >>%paste

#===================================================================
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

#===================================================================

cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Speleothem_Development/PRYSM/Tester_Files/

ST=np.load('HIDDEN_TEMP.npy')
SP=np.load('HIDDEN_PRCP.npy')
Sd18Op=np.load('HIDDEN_d18OP.npy')
Sd18Os=np.load('HIDDEN_d18OS.npy')


H1=np.load('HIDDEN_PDRIP_T_025.npy')
H2=np.load('HIDDEN_PDRIP_T_05.npy')
H3=np.load('HIDDEN_PDRIP_T_1.npy')
H4=np.load('HIDDEN_PDRIP_T_5.npy')


H1=H1-np.mean(H1)
H2=H2-np.mean(H2)
H3=H3-np.mean(H3)
H4=H4-np.mean(H4)

H1f, H1psd_mt, H1nu = tsa.multi_taper_psd(H1, Fs=1.0,adaptive=False, jackknife=False)
H2f, H2psd_mt, H2nu = tsa.multi_taper_psd(H2, Fs=1.0,adaptive=False, jackknife=False)
H3f, H3psd_mt, H3nu = tsa.multi_taper_psd(H3, Fs=1.0,adaptive=False, jackknife=False)
H4f, H4psd_mt, H4nu = tsa.multi_taper_psd(H4, Fs=1.0,adaptive=False, jackknife=False)


ST=ST-np.mean(ST)
SP=SP-np.mean(SP)
Sd18Op=Sd18Op-np.mean(Sd18Op)
Sd18Os=Sd18Os-np.mean(Sd18Os)

STf, STpsd_mt, STnu = tsa.multi_taper_psd(ST, Fs=1.0,adaptive=False, jackknife=False)
SPf, SPpsd_mt, SPnu = tsa.multi_taper_psd(SP, Fs=1.0,adaptive=False, jackknife=False)
Sd18Opf, Sd18Oppsd_mt, Sd18Opnu = tsa.multi_taper_psd(Sd18Op, Fs=1.0,adaptive=False, jackknife=False)
Sd18Osf, Sd18Ospsd_mt, Sd18Osnu = tsa.multi_taper_psd(Sd18Os, Fs=1.0,adaptive=False, jackknife=False)

t=np.arange(1000,2005,1)

#===================================================================

# Calculate Spectral Slope for each cave.

ifit = np.where((H1f > 1./200.) & (H1f< 1./2.))
slope1 = np.polyfit(log(H1f[ifit]),log(H1psd_mt[ifit]),1)
slope2 = np.polyfit(log(H2f[ifit]),log(H2psd_mt[ifit]),1)
slope3 = np.polyfit(log(H3f[ifit]),log(H3psd_mt[ifit]),1)
slope4 = np.polyfit(log(H4f[ifit]),log(H4psd_mt[ifit]),1)

# 500:
# In [11]: slope1
# Out[11]: array([ 0.04648435,  2.33212079])
# 
# In [12]: slope2
# Out[12]: array([-0.09851548,  0.75640274])
# In [14]: slope3
# Out[14]: array([-0.45501324, -1.22535258])
# In [17]: slope4
# Out[17]: array([-0.84986353, -4.25305739])
# 
# #200:
# In [25]: slope1
# Out[25]: array([ 0.04087977,  2.32402961])
# 
# In [26]: slope2
# Out[26]: array([-0.1141673 ,  0.73376226])
# 
# In [27]: slope3
# Out[27]: array([-0.4919703 , -1.27883897])
# 
# In [28]: slope4
# Out[28]: array([-0.88200078, -4.29960644])


#======================================================================================


plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.rcParams['font.family']='serif'

ax=plt.figure()
plt.plot(t,H1,label=r'$\tau=0.25',color='Navy')
plt.plot(t,H2,label=r'$\tau=0.5',color='MediumBlue')
plt.plot(t,H3,label=r'$\tau=1',color='SlateBlue')
plt.plot(t,H4,label=r'$\tau=5',color='Purple')

plt.title(r'Hidden Cave, NM: Effects of Groundwater Transit Time, $\tau$',fontsize=16)
plt.xlabel(r'Years', fontsize=11)
plt.ylabel(r'$\delta^{18}O_C$', fontsize=11)
#plt.spines["top"].set_visible(False)  
#ax.spines["right"].set_visible(False)  

plt.legend(loc=3,fontsize=11,frameon=False)

#======================================================================================

# Plot Spectra for all
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True
plt.rcParams['font.family']='serif'

ax1=plt.subplot(2,1,1)
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False) 
ax1.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.loglog(STf,STpsd_mt, label='T',color='red')
plt.loglog(SPf,SPpsd_mt, label='P',color='blue')
plt.loglog(Sd18Opf,Sd18Oppsd_mt, label=r'$\delta^{18}O_P$',color='DodgerBlue')
plt.xlabel(r'Period (Years)')
plt.ylabel(r'PSD')
plt.title(r'ENVIRONMENT',fontsize=11, color='gray')
ax1.minorticks_off()

# Change the x-tick labels to something reasonable
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.grid('on',axis='y',color='DimGray')
plt.xlim([1./550.,0.4])
plt.ylim([1e-3, 1e1])
plt.legend(loc=3,ncol=3,fontsize=14,frameon=False)
#======================================================================================


ax2=plt.subplot(2,1,2)

ax2.spines["top"].set_visible(False)  
ax2.spines["right"].set_visible(False) 
ax2.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
plt.loglog(H1f,H1psd_mt,label=r'$\tau=0.25',color='Navy',linewidth=2)
plt.loglog(H2f,H2psd_mt,label=r'$\tau=0.5',color='OrangeRed',linewidth=2)
plt.loglog(H3f,H3psd_mt,label=r'$\tau=1',color='LightSeaGreen',linewidth=2)
plt.loglog(H4f,H4psd_mt,label=r'$\tau=5',color='DarkMagenta',linewidth=2)
plt.tick_params(axis="both", which="both", bottom="on", top="off",  
                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
pertick=np.array([500, 200, 100 ,50, 20,10,8, 6,4, 2])
pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

xtick=(1.0/pertick)

plt.xticks(xtick, pertick_labels, size='small')
plt.minorticks_off()
plt.grid('on',axis='y',color='DimGray')
plt.title(r'SENSOR SIGNAL',fontsize=11, color='gray')
plt.xlim([1./550.,0.4])
#plt.ylim([1e-4, 1])
plt.xlabel(r'Period (Years)', fontsize=11)
plt.ylabel(r'PSD', fontsize=11)
#plt.spines["top"].set_visible(False)  
#ax.spines["right"].set_visible(False)  

plt.legend(loc=3,fontsize=14,frameon=False)
plt.suptitle(r'Effects of Groundwater Transit Time, $\tau$ (Hidden Cave, NM)',fontsize=16)

plt.show()
