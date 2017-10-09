# Sylvia Dee <sdee@usc.edu>
# OBSERVATION SUB-MODEL
# Function BAM 'psm_obs_banded'
# Modified [sdee] 10/28/2013

#====================================================================

# BAM_SIMUL_PERTURB for Python
# Original Code in Matlab by Maud Comboul <comboul@usc.edu>
# Revised for Python by <sdee@usc.edu>
# 6/12/14

'''
Observation Model: BAM

Use: generate an ensemble of age perturbed data.

See www.clim-past-discuss.net/9/6077/2013/ for a detailed description of the model.
The time series in X are automatically flipped to range from most recent to oldest measurements when the intput t is given in increasing order.

[Xp tp tmc] = BAM_simul_perturb(X,t) will generate an ensemble of 1000 age models randomly following
a Poisson process with rate parameter theta=0.05 used to perturb data X

[Xp tp tmc] = BAM_simul_perturb(X,t,**kwargs) will perturb data X  with the model specified in
the model structure.


INPUT
    X: data (vector or matrix n*p)
    t: chronology for data X (n*1)

KEYWORD-ARGUMENTS:
    ns: number of samples
    name: 'poisson' or 'bernoulli'
    param: probability of growth band being perturbed (default: prob of missing band = prob of doubly-counted band = 0.05)
         if param is a single argument, then the perturbations are symmetric (prob of missing band = prob of doubly-counted band)
         if param = [a1 a2] and a1 neq a2 the model is asymmetric
                        a1 = prob(missing layer)
                        a2 = prob(layer counted multiple times)
    resize: do not resize: 0 (default), resize to shortest sample: -1, resize to longest sample: 1
    tm: if a time model is provided, the code returns the corresponding perturbed data

OUTPUT
    Xp: realizations of age-perturbed data matrix of size tn*p*ns (could be 2 or 3d)
    tp: new chronology tn*1
    tmc: corresponding ensemble of time-correction matrices (tn*p*ns) to map realizations in Xp back to the original data X (2=insert nan, 0=remove double band) (2 or 3d)
    where tn is the chronology length = n (default), shortest sample or longest sample
    depending on the chosen resizing option.

    '''

#======================================================================
#  0.0 Initialization
#======================================================================

import numpy as np
import random

#======================================================================
# 1.0 Function Bam Simul Perturb
#====================================================================

def bam_simul_perturb(X,t,param=[0.01,0.01],name='poisson',ns=1000,resize=0):


# 1.0 Prep Data
#
# Use only for multi-column data.

#   if np.size(X,0) < np.size(X,1):
#       X = X.T

    t=t[:]
    if np.mean(np.diff(t))>0:
            isflipped = 1
            X = np.flipud(X)
            t = np.flipud(t)
    else:
            isflipped = 0

            X=X.T
            X=X.T

    n,p=np.shape(X)

    if len(param)==1:
        param[1] = param[0]


#====================================================================
# 2.0 Generate an ensemble of time perturbation models
#====================================================================

    if not 'tm' in locals():
        tm = np.ones((n,p,ns))
        if name.lower()=='poisson'.lower():  # compares strings for equality ignoring case
            # apply poisson model
            for nn in range(ns):

                num_event_mis = np.random.poisson(param[0]*n,size=(p,1)) # poisson model for missing bands
                num_event_dbl = np.random.poisson(param[1]*n,size=(p,1)) # poisson model for layers counted multiple times

                        #print(num_event_mis[0][0])
                for ii in range(p):
                    #nnew=n-1
                    jumps = np.random.choice(n-1,num_event_mis[ii][0])+1            # place events uniformly on {2,...,n}
                    tm[jumps,ii,nn] = tm[jumps,ii,nn]-1                 # remove 1 at jump locations
                    jumps = np.random.choice(n-1,num_event_dbl[ii][0])+1
                    tm[jumps,ii,nn] = tm[jumps,ii,nn]+1                 # add 1 at jump locations


        elif name.lower()=='bernoulli'.lower():
            for nn in range(ns):

                # apply bernoulli model
                #numpy.random.binomial(n, p, size=None)

                tm=tm-np.random.binomial(1,param[0],size=(n,p))
                tm=tm+np.random.binomial(1,param[1],size=(n,p))

                #tm = tm-binornd(1,param[0],[n,p]) # Binomial model for missing bands
                #tm = tm+binornd(1,param[1],[n,p]) # Binomial model for layers counted multiple times


        else:
            print('Unknown age model ; acceptable inputs are ''poisson'' and ''bernoulli''')

#====================================================================
# generate age perturbed data Xp
# expand length of Xp and tXp if resizing is required
#====================================================================

    if resize==1:
        t_ext = np.ceil(2*param(2)*n)
        tn = n + t_ext
        X = np.concatenate(X, np.nan(t_ext,p))

        dt = t(2)-t(1)
        time_ext = t[-1]+np.arange(dt,t[-1],dt)+t_ext*dt.T      # END--not python.
        tp = np.concatenate(t, time_ext.T)

    else:
        tn=n
        tp=t

#====================================================================
#====================================================================

    Xp = np.empty([tn,p,ns])
    for i in range(tn):
        for j in range(p):
            for k in range(ns):
                Xp[i,j,k]=np.nan



    Tmax = 0
    Tmin = n

    tmc=np.ones((tn,p,ns))
    for nn in range(ns):

        for ii in range(p):
            xcount=0
            Xcount=0
            tt=0
            while tt < n:

                if tm[tt,ii,nn]==0:   # remove band
                    Xcount=min(Xcount+1,tn)
                    tmc[xcount-1,ii,nn]=tmc[xcount-1,ii,nn]+1
                elif   tm[tt,ii,nn]==2:   # insert double band
                    Xp[xcount-1,ii,nn]=X[Xcount-1,ii]
                    tmc[xcount-1,ii,nn]=tmc[xcount-1,ii,nn]-1
                    xcount=min(tn,xcount+1)
                #print(Xcount,xcount)

                Xp[xcount-1,ii,nn]=X[Xcount-1,ii]
                xcount=min(tn,xcount+1)
                Xcount=min(tn,Xcount+1)
                tt=tt+1
#np.where(~numpy.isnan(Xp)[0][-1]
            k=np.where(~np.isnan(Xp[:,ii,nn]))[0][-1] #gives last element
            if k > Tmax:
                Tmax=k

            if k < Tmin:
                Tmin=k

#====================================================================
# Resize all final data
#====================================================================

    # crop output size to shortest sequence
    if resize==-1:
        Xp = Xp[1:Tmin,:,:]
        tp=tp[1:Tmin]


    # expand output size to longest (non-nan) sequence
    if resize==1:
        Xp = Xp[1:Tmax,:,:]
        tp=tp[1:Tmax]


    if X.shape[1]==1:
        Xp = np.reshape(Xp,(n,ns))
        tmc = np.reshape(tmc,(n,ns))


    if isflipped == 1:
        if X.shape[1]==1:
            Xp = np.flipud(Xp)
        else:
            for nn in range(ns):
                Xp[:,:,nn] = np.flipud(Xp[:,:,nn])
                tmc[:,:,nn] = np.flipud(tmc[:,:,nn])


        tp = np.flipud(tp)

#====================================================================
# Return simulated age models
#====================================================================

    return tp, Xp, tmc
