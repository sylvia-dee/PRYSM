
# Sylvia Dee <sylvia_dee@brown.edu>
# PRYSM
# PSM for Lacustrine Sedimentary Archives
# ARCHIVE MODEL: Bioturbation of Sediments
# Function 'bioturbation'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
# Modified 04/2/2018 <sylvia@ig.utexas.edu>

#====================================================================

def bioturbation(abu,iso,mxl,numb):
    '''
    function [oriabu,bioabu,oriiso,bioiso] = turbo2(abu,iso,mxl,numb)
    # The MATLAB program TURBO2 can be used to simulate the effects of
    # bioturbation on single sediments particles.
    #
    # Converted to PYTHON by Sylvia Dee, <sylvia_dee@brown.edu>
    #
    # Working Variables:
    # ABU = series of abundances of carrier type 1 down core
    # ISO = series of isotopes measured on carrier
    # MXL = series of mixed layer thicknesses down core
    # NUMB = number of carriers to be measured

    # Outputs:
    # ORIABU = original abundances of both carrier types 1 and 2
    # BIOABU = bioturbated abundances of both carriers types 1 and 2
    # ORIISO = original isotope signature of both carrier  types 1 and 2
    # BIOISO = bioturbated isotope signature of both carrier types 1 and 2
    
    # The first column contains the age of the sediment. Note that the array 
    # starts with the oldest layer, which is deposited first. In our example, 
    # the age (in calendar kyr BP) decreases linearly with time, i.e. by 0.2 kyr
    # every centimeter corresponding to a sedimentation rate of 5 cm kyr−1. The 
    # second column contains the mixed layer thickness in centimeters; in the 
    # example the mixed layer thickness is constant at 10 cm. The third column 
    # contains the abundance of Species 1. The fourth column contains the isotope 
    # signal from Species 1 (in our example oxygen isotope values are in ‰ vs. PDB).

    # Martin Trauth 18 July 2012
    # Adapted for Python/PRYSM Paloclimate modeling toolbox March 2017
    '''
 
    import numpy as np
    import scipy
    from scipy import integrate
    from scipy import stats
    import matplotlib.pyplot as plt
    from math import pi, sqrt, exp

    nrows = int(np.max(mxl)+0)
    ncols = int(np.max(abu)+50)
    sedabu = np.zeros((nrows,ncols))
    sedabu.fill(np.nan)
    sediso = np.zeros((nrows,ncols))
    sediso.fill(np.nan)
    #new=np.zeros((len(abu)+10,ncols))
    #h = waitbar(0,'Mixing Process ...');
    mxl=np.array(mxl,dtype=int)

    ni=np.zeros((sedabu.shape))
    ni.fill(np.nan)
    ns=np.zeros((sedabu.shape))
    ns.fill(np.nan)


    for i in range(len(abu)):
        #waitbar(i/length(abu))
        rncols=np.random.permutation(ncols)
        #sedabu[sedabu.shape[0]+1,0:ncols] = 2*np.ones((1,ncols))
        sedabu=np.append(sedabu,np.ones((1,ncols))*2,axis=0) 
        #sedabu(size(sedabu,1),1:abu(i)) = ones(1,abu(i));
        sedabu[len(sedabu[:,0])-1,0:int(abu[i])]=np.ones((1,abu[i]))
        #sedabu(size(sedabu,1),:) = sedabu(size(sedabu,1),rncols);
        sedabu[len(sedabu[:,0])-1,:] = sedabu[len(sedabu[:,0])-1,rncols]
        # sedabu(size(sedabu,1)+1,1:ncols) = 2*ones(1,ncols);
        #sediso[len(sediso[:,0])-1,0:ncols] = iso[i]*np.ones((1,ncols))
        #sediso(size(sediso,1)+1,1:ncols) = iso(i)*ones(1,ncols);
        sediso=np.append(sediso,iso[i]*np.ones((1,ncols)),axis=0)
        for j in range(ncols):
            z = np.random.permutation(int(mxl[i]))
            ns[0:mxl[i],j] = sedabu[(len(sedabu[:,0])-mxl[i]):len(sedabu[:,0]),j]
            sedabu[(len(sedabu[:,0])-mxl[i]):len(sedabu[:,0]),j] = ns[z,j]
            ni[0:mxl[i],j] = sediso[(len(sediso[:,0])-mxl[i]):len(sediso[:,0]),j]
            sediso[(len(sediso[:,0])-mxl[i]):len(sediso[:,0]),j] = ni[z,j]

    #      sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),j) = ni(z,j);

    # ======================================================================

    # clear ns ni i j mxl z (SDEE IS THIS IMPORANT?)

    # initialize other arrays:
    oriabu=np.zeros((len(abu),2))
    bioabu=np.zeros((len(abu),2))
    oriiso=np.zeros((len(abu),2))

    # fill final outputs
    sedabu = sedabu[nrows::,:]
    sediso = sediso[nrows::,:]

    oriabu[:,0] = abu
    oriabu[:,1] = ncols-abu

    x1=np.where(sedabu==1.)
    x2=np.where(sedabu==2.)
    sum_array1=np.zeros((sedabu.shape))
    sum_array2=np.zeros((sedabu.shape))
    sum_array1[x1]=1.
    sum_array2[x2]=1.

    bioabu[:,0] = np.sum(sum_array1,axis=1)
    bioabu[:,1] = np.sum(sum_array2,axis=1)
    oriiso[:,0] = iso
    oriiso[:,1] = iso

    bioiso1 = sediso
    bioiso2 = sediso

    # ======================================================================

    #for i in range(len(sedabu[:,0]))
    # x=np.where((sedabu!=1.))
    # y=np.where((sedabu!=2.))

    #bioiso1[sedabu!=1.]=np.nan
    #bioiso2[sedabu!=2.]=np.nan

    def nanfunc1(a,b):
    	return np.where(b!=1., np.nan, a)


    test1=nanfunc1(bioiso1,sedabu)

    def nanfunc2(a,b):
    	return np.where(b!=2., np.nan, a)

    test2=nanfunc2(bioiso2,sedabu)

    bioiso1=test1
    bioiso2=test2

    # for i in range(239):
    # 	for j in range(540):
    # 		if (sedabu[i,j]!=1.):
    # 			bioiso1[i,j]=np.nan
    # 		elif (sedabu[i,j]!=2.):
    # 			bioiso2[i,j]=np.nan


    # ======================================================================

    # bioiso1[x].shape
    # bioiso1[x] = np.nan
    # bioiso2[y] = np.nan

    # ======================================================================


    biopart1 = np.zeros((bioiso1.shape))
    biopart1.fill(np.nan)
    biopart2 = np.zeros((bioiso2.shape))
    biopart2.fill(np.nan)


    # for i in range(len(abu)):
    #     temp=np.where(np.isnan(bioiso1[i,:]))
    #     bioiso1[i,temp]=0
    #     biopart1[i,0:int(bioabu[i,0])] = bioiso1[i,np.sum(bioiso1[i,temp])]
    #     temp2=np.where(np.isnan(bioiso2[i,:]))
    #     bioiso2[i,temp2]=0
    #     biopart2[i,0:int(bioabu[i,1])] = bioiso2[i,np.sum(bioiso2[i,temp2])]


    for i in range(len(abu)):
    	x=np.where(~np.isnan(bioiso1[i,:]))
    	biopart1[i,0:int(bioabu[i,0])] = bioiso1[i,x]
    	y=np.where(~np.isnan(bioiso2[i,:]))
    	biopart2[i,0:int(bioabu[i,1])] = bioiso2[i,y]



    biopart1 = biopart1[:,0:numb]
    biopart2 = biopart2[:,0:numb]

    bioiso=np.zeros((len(abu),2))
    for i in range(len(abu)):
        bioiso[i,0] = np.nanmean(biopart1[i,:])
        bioiso[i,1] = np.nanmean(biopart2[i,:])


    return oriabu,bioabu,oriiso,bioiso


# =====================================================================
# mxltext = str(np.mean(mxl))
# numbtxt = str(numb)

# # PLOT
# plt.style.use('ggplot')

# fig2 = plt.figure(1)

# fig2.add_subplot(1,2,1)
# plt.plot(years,oriiso[:,0],'deepskyblue',label=r'Original TEX86',linewidth=3.0)
# plt.plot(years,bioiso[:,0],'peru',label=r'BioCarrier 1')
# plt.plot(years,bioiso[:,1],'firebrick',label=r'BioCarrier 2')
# plt.grid(True)
# plt.gca().invert_yaxis()

# #plt.title(r'Isotopes of Carriers 1+2,'+  mxltext + 
# #    ' cm Mixed Layer, '+ numbtxt +' Carriers')
# plt.legend()

# fig2.add_subplot(1,2,2)
# plt.fill_between(years,TEX_low,TEX_high,color='gray',label=r'TEX86 Calibration Error')
# plt.plot(years,oriiso[:,0],'deepskyblue',label=r'Original TEX86',linewidth=3.0)
# plt.plot(years,bioiso[:,0],'peru',label=r'BioCarrier 1')
# plt.plot(years,bioiso[:,1],'firebrick',label=r'BioCarrier 2')
# plt.grid(True)
# plt.gca().invert_yaxis()

# plt.suptitle(r'$TEX_{86}$ in Carriers 1+2,'+  mxltext + 
#     ' cm Mixed Layer, '+ numbtxt +' Carriers')
# plt.legend()
# plt.show()

# ======================================================================

