# Sylvia Dee <sdee@usc.edu>
# PRYSM
# Post-Processing and Figures
# Figure 6: Age Uncertainties
# Script to reproduce paper results
# Modified [sdee] 01/30/2015
#===================================================================

# NOTE TO USER: THIS SCRIPT IS DESIGNED TO BE RUN IN IPYTHON CONSOLE
# highlight, copy all text, >>%paste

#===================================================================
# Initialization

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

#======================================================================================
# Load Climate and Isotope Variables

cd /home/geovault-02/sdee/Forward_Modeling_Toolbox/Cellulose_Development/FigureVars/
TT=np.load('LaSelva_T.npy')
TP=np.load('LaSelva_P.npy')
Td18Op=np.load('LaSelva_d18OP.npy')
Td18Os=np.load('LaSelva_d18OS.npy')
TV=np.load('LaSelva_d18OV.npy')
TRH=np.load('LaSelva_RH.npy')
TC=np.load('LaSelva_d18Ocell.npy')

# set time axis
t=np.arange(1000,2005,1)

#======================================================================================

# Load Age Uncertainties

# Run Bam_SIMUL
#X=TC.reshape(1005,1)
#tp, Xp, tmc =bam_simul_perturb(X,t,param=[0.02,0.02])
Xp=np.load('LaSelva_AgePerturbed_Network.npy')

# Plot some BAM with 4% total error
from scipy.stats.mstats import mquantiles
q1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
q2=t

#======================================================================================
figure(1)
plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations, CI',facecolor='gray',alpha=0.5)
plt.plot(t,TC,color='DarkGreen')
p2 = Rectangle((0, 0), 1, 1, fc="gray")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="DarkGreen")
legend([p1,p2],[r'$\delta^{18}O_{C}$','1000 Age-Perturbed Realizations, [2.5 97.5 CI]'],loc=3,fontsize=14,frameon=False)
plt.title(r'La Selva, Costa Rica: \texttt{BAM} simulated Age Perturbations [4$\%$ Error Rate]')
plt.grid('on')
plt.xlim([1000,2010])
plt.ylim([26,28.5])
plt.xlabel(r'Year')
plt.ylabel(r'$\delta^{18}O_{C}$')
plt.show()

#======================================================================================
# Correlation between La Selva and NINO34
#======================================================================================

nino=np.load('T_NINO34.npy')

# Compute anomalies
nino_anom=nino-np.mean(nino)
#Xp
Xpm=Xp-np.mean(Xp, axis=0)
TC=TC-np.mean(TC)

from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(nino_anom,TC)

BAMslopes=np.zeros((1000))
BAMints=np.zeros((1000))
BAMr=np.zeros((1000))
ps=np.zeros((1000))
stds=np.zeros((1000))
for i in range(1000):
    BAMslopes[i],BAMints[i],BAMr[i],ps[i],stds[i]= stats.linregress(nino_anom,Xpm[:,i])


#======================================================================================
# Plot Values
x=np.linspace(-3,3)

from scipy.stats.mstats import mquantiles
slopes=mquantiles(BAMslopes,prob=[0.025,0.975],axis=0)
ints=mquantiles(BAMints,prob=[0.025,0.975],axis=0)

line_0=slope*x+intercept
line_low=slopes[0]*x+ints[0]
line_high=slopes[1]*x+ints[1]

figure(2)
plt.plot(x,line_0,color='Blue',linestyle='-',label='Regression: Unperturbed Time Series (R= - 0.52)')
plt.fill_between(x,line_low,line_high,facecolor='Gray')
plt.scatter(nino_anom,TC,color='ForestGreen',marker='v',label=r'$\delta^{18}O_{C}$ vs. NINO34')
plt.xlim([-3,3])
plt.title('La Selva Age Perturbed Cellulose Data: Regression against NINO34 Index, $4\%$ Dating Uncertainty',y=1.05)
plt.xlabel('NINO3.4 SST Anomaly [K]')
plt.ylabel(r'La Selva Tree Ring Cellulose $\delta^{18}O_{C}$ Anomaly')
plt.grid('on')
p3 = Rectangle((0, 0), 1, 1, fc="gray")
p2 = Rectangle((0, 0), 0.7, 0.7, fc="Blue")
p1 = Rectangle((0, 0), 0.5, 0.5, fc="ForestGreen")
legend([p3,p2,p1],[r'Regression for Age-Perturbed Realizations, [2.5 97.5 CI]','Unperturbed Regression ($R=-0.52, p<<0.01$)','$\delta^{18}O_{C}$ vs. NINO34 ST',],loc=3,fontsize=14,frameon=False)
#plt.legend(loc=3,fontsize=14,frameon=False)
plt.show()
#======================================================================================