# ***************************************

# Sylvia Dee <sylvia_dee@brown.edu>
# Script to Read Input Files and Run Leaf Wax Model

#======================================================================
# LEAF WAX PSM: Bronwen Konecky (original in NCL)
#               Modified 7/21/17 Sylvia Dee
#======================================================================

# PSEUDO CODE

# 1. read in dD precip, soil/xylem, PFT
# 2. regrid and apply epsilon values to grid based on PFT values (E29 values) (14 x lat x lon)
# 3. calculate percentage PFT f_(plant type), Etot ~ weighted average eps. value for grid cell
# 4. apply to dD wax

# OR Input Etot as netcdf and apply to CAM, regrid for SPEEDY/IsoGSM

#======================================================================
# 0.0 Initialization
#======================================================================
from pylab import *
import sys, os, cdtime, cdutil, cdms2, vcs, MV2, time, datetime
import numpy as np
import matplotlib.pyplot as plt
from math import exp

#======================================================================
# 1.0 Load Data for each Climate Model
#======================================================================

# start_time=cdtime.comptime(1980)
# print 'start_time = ', start_time
# end_time=cdtime.comptime(2005)
# print 'end_time =', end_time

#======================================================================
# CAM5 (Load Data)
#======================================================================

dir='/disk/staff/sylvia/LEAFWAX/'

file_name=dir + 'f.e13.FiPIPDC5CN.f19_f19.ctrl.cam.h0.1850-2014.prec.nc'

#append file name

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

# MANUALLY PULL DATES

CAM5_PR = f('P',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))
#cdutil.setTimeBoundsMonthly(CAM5_PR)
CAM5_PR=CAM5_PR[1560:-120,:,:]
#======================================================================

CAM5_dD_pr = f('dD',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))
CAM5_dD_pr=CAM5_dD_pr[1560:-120,:,:]

#======================================================================
times=np.arange(1560.0,1860.0,1)
newTimeAxis = cdms2.createAxis(times, id='time')
newTimeAxis.units = 'months since 1850'
newTimeAxis.designateTime()
CAM5_PR.setAxis(0,newTimeAxis)

cdutil.setTimeBoundsMonthly(CAM5_PR)

#======================================================================
times=np.arange(1560.0,1860.0,1)
newTimeAxis = cdms2.createAxis(times, id='time')
newTimeAxis.units = 'months since 1850'
newTimeAxis.designateTime()
CAM5_dD_pr.setAxis(0,newTimeAxis)

cdutil.setTimeBoundsMonthly(CAM5_dD_pr)
#======================================================================
file_name=dir + 'f.e13.FiPIPDC5CN.f19_f19.ctrl.clm2.h0.dDxylem.1850-2014.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

CAM5_dD_xy = f('dDxylem',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))
#cdutil.setTimeBoundsMonthly(CAM5_dD_xy)
CAM5_dD_xy=CAM5_dD_xy[1560:-120,:,:]
#======================================================================

times=np.arange(1560.0,1860.0,1)
newTimeAxis = cdms2.createAxis(times, id='time')
newTimeAxis.units = 'months since 1850'
newTimeAxis.designateTime()
CAM5_dD_xy.setAxis(0,newTimeAxis)

cdutil.setTimeBoundsMonthly(CAM5_dD_xy)

#======================================================================

# CAM5 dD of Precip
CAM_dD_pr_ANN=cdutil.YEAR(CAM5_dD_pr)

# CAM5 dD of Xylem Water
CAM_dD_xy_ANN=cdutil.YEAR(CAM5_dD_xy)

# CAM5 Precip Amount
CAM_pr_ANN=cdutil.YEAR(CAM5_PR)

#======================================================================
# EPSILON VALUES
#======================================================================


#======================================================================
# 2. ETOT ~ read for replication/sanity check
#======================================================================

dir = '/disk/staff/sylvia/LEAFWAX/'

file_name=dir + 'Etot.year.2000.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#file_name=dir +'Etot.year.1850.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#append file name

#file_name=dir+'Etot.year.-19000.trace21k.pft.lgm.nc'

# LGM: the land mask changes between the two time periods, so a difference map of predicted 
# dDc29 LGM-2000 will probably give you errors. For this figure, I'd recommend just masking out 
# anything that has a missing value either during the LGM or during the year 2000. That will take c
# are of points like Sundaland that have data during the LGM but not 2000, and points like under 
# the Laurentide ice sheet that have data during 2000 but not during the LGM. 


f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

ETOT = f('Etot',longitude=(0,360), latitude = (-90., 90.))
# ETOT_1850= f('Etot',longitude=(0,360), latitude = (-90., 90.))

# create mask where Etot=0
ETOT_mask=np.ma.masked_equal(ETOT, -999.)
# SD seems to already be masked coming out of NCL.

#======================================================================
# 3. PFT VEGETATION
#======================================================================

dir = '/disk/staff/sylvia/LEAFWAX/'

file_name=dir + 'surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'


f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

PCT = f('PCT_PFT',longitude=(0,360), latitude = (-90., 90.))

PCT_avg=np.mean(PCT,axis=0)

PCT_mask=np.ma.masked_equal(ETOT, -999.)

#======================================================================
# SET THESE VALUES IF YOU WANT THE RAW DATA, EPSILON HAS NO UNCERTAINTY

E29_sh = -99.
E29_tr = -121.
E29_fo = -128.
E29_c3 = -149.
E29_c4 = -132.

#======================================================================

# THIS IS FOR MCA, ASSUMING EPSILON ESTIMATES ARE UNCERTAIN

SH_MCA=np.linspace(-67, -130,1000)
TR_MCA=np.linspace(-99,	-143,1000)
FO_MCA=np.linspace(-98,	-158,1000)
C3_MCA=np.linspace(-120,-177,1000)
C4_MCA=np.linspace(-107,-157,1000)

from numpy import random

# INITIALIZE MCA ARRAY

Etot_PCT=np.zeros((1000,96,144))

# COMPUTE NEW ETOT VALUE 1000 TIMES, SAMPLING UNCERTAINTIES IN EPSILON:

for k in range(1000):

	E29_sh = np.random.choice(SH_MCA,1)
	E29_tr = np.random.choice(TR_MCA,1)
	E29_fo = np.random.choice(FO_MCA,1)
	E29_c3 = np.random.choice(C3_MCA,1)
	E29_c4 = np.random.choice(C4_MCA,1)

	E29_full=np.zeros((17))

	E29_full[0]  = -999.
	E29_full[1]  = E29_tr
	E29_full[2]  = E29_tr
	E29_full[3]  = E29_tr
	E29_full[4]  = E29_tr
	E29_full[5]  = E29_tr
	E29_full[6]  = E29_tr
	E29_full[7]  = E29_tr
	E29_full[8]  = E29_tr
	E29_full[9]  = E29_sh
	E29_full[10] = E29_sh
	E29_full[11] = E29_sh
	E29_full[12] = E29_c3
	E29_full[13] = E29_c3
	E29_full[14] = E29_c4
	E29_full[15]  = -999.
	E29_full[16] = -999.


	# ; Needleleaf evergreen tree - temperate
	# ; Needleleaf evergreen tree - boreal
	# ; Needleleaf deciduous tree - boreal
	# ; Broadleaf evergreen tree - tropical
	# ; Broadleaf evergreen tree - temperate
	# ; Broadleaf deciduous tree - tropical
	# ; Broadleaf deciduous tree - temperate
	# ; Broadleaf deciduous tree - boreal
	# ; Broadleaf evergreen shrub - temperate
	# ; Broadleaf deciduous shrub - temperate
	# ; Broadleaf deciduous shrub - boreal
	# ; C3 arctic grass
	# ; C3 grass
	# ; C4 grass

	E29gl=PCT_avg


	E29gl[0,:,:]=-999.
	E29gl[15,:,:]=-999.
	E29gl[16,:,:]=-999.

	E29gl=np.ma.masked_equal(E29gl, -999.)
	E29_full=np.ma.masked_equal(E29_full, -999.)

	veg_pct=E29gl
	veg_sum=np.sum(veg_pct,axis=0)

	for i in range(17):
		veg_pct[i,:,:]= veg_pct[i,:,:] / veg_sum[:,:]


	for i in range(17):
	    E29gl[i,:,:]= veg_pct[i,:,:] * E29_full[i]


	Etot_PCT[k,:,:]=np.sum(E29gl,axis=0)


#======================================================================
# SET UP CALCULATION, PICK VARIABLE 
#======================================================================
cd /disk/staff/sylvia/LEAFWAX/
Etot_PCT=np.load('ETOT_PCT_MCA.npy')
# extract percentiles of Etot

Q1=np.percentile(Etot_PCT,2.5,axis=0)
Q2=np.percentile(Etot_PCT,97.5,axis=0)
#======================================================================
# SET TO MAP AVERAGE

CAM_5_mask=np.ma.masked_greater(CAM5_dD_pr,1e+9)
CAM5_mask=CAM_5_mask.reshape(25,12,96,145)
CAM5_ANN_pr=np.mean(CAM5_mask,axis=(0,1))

CAM_5_mask=np.ma.masked_greater(CAM5_dD_xy,1e+9)
CAM5_mask=CAM_5_mask.reshape(25,12,96,145)
CAM5_ANN_xy=np.mean(CAM5_mask,axis=(0,1))

#======================================================================

# RESAMPLING/MONTE CARLO PROCESS:
# careful of mask

# mask=ma.getmaskarray(ETOT_mask)
# ma.MaskedArray.get_fill_value(ETOT_mask)

# index=np.where(mask==False)

# # create array of ETOTs from random distribution
# from numpy import random

# epsilon_mca=np.linspace(-140,-90,1000)
# ETOT_MCA=np.ma.zeros((1000,96,145))

# for i in range(1000):
#     ETOT_MCA[i][index]=np.random.choice(epsilon_mca,1)

# for i in range(1000):
# 	ETOT_MCA[i,:,:]=np.ma.array(ETOT_MCA[i,:,:],mask=mask)

#======================================================================

# clip extra longitude point
CAM5_ANN_pr=CAM5_ANN_pr[:,0:144]
CAM5_ANN_xy=CAM5_ANN_xy[:,0:144]

#======================================================================

# Calculate dD wax 1000 times for each grid cell, then get the distribution.

CAM_WAX_pr=np.zeros((1000,96,144))
for i in range(1000):
	CAM_WAX_pr[i,:,:]=((CAM5_ANN_pr + 1000.) * (Etot_PCT[i,:,:]/1000. + 1.)) - 1000.


# mask unrealistic high values, artifact of shortcomings of numpy.zeros

CAM_WAX_MCA=np.ma.masked_greater(CAM_WAX_pr,10.)

# calculate standard deviation of dDwax in each grid cell.

CAM_WAX_STD=np.std(CAM_WAX_MCA,axis=0)

#======================================================================

# Use Quantiles instead:

CAM_WAX_pr_25=((CAM5_ANN_pr + 1000.) * (Q1/1000. + 1.)) - 1000.
CAM_WAX_pr_97=((CAM5_ANN_pr + 1000.) * (Q2/1000. + 1.)) - 1000.

# 95% CI difference in dD calculated:
CAM_WAX_pr=CAM_WAX_pr_97-CAM_WAX_pr_25

CAM_WAX_MCA=CAM_WAX_pr

CAM_WAX_STD=CAM_WAX_MCA

#======================================================================
# MAPS // PLOTTING
#======================================================================

C_lat=np.array(CAM5_PR.getLatitude())
C_lon=np.array(CAM5_PR.getLongitude())
C_lon=C_lon[0:144]


# READ IN LAND MASK
dir = '/disk/staff/sylvia/LEAFWAX/'
file_name=dir + 'landmask.f.e13.FiPIPDC5CN.f19_f19.ctrl.cam.h0.1850-2014.prec.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

landmask = f('landmask',longitude=(0,359), latitude = (-90., 90.))

thresh=0.5
land_threshold_indices = landmask > thresh
ocean_threshold_indices = landmask < thresh

landmask[land_threshold_indices] = 1.0
landmask[ocean_threshold_indices] = 0.0

CAM_WAX_land=CAM_WAX_STD*landmask
CAM_MAP_mask=np.ma.masked_inside(CAM_WAX_land, -0.1, 0.1)

#======================================================================

from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

CAM_MAP2,C_lon2 = shiftgrid(180.,CAM_MAP_mask,C_lon,start=False)

vmin=25.
vmax=55.
levels=np.arange(25,55,1)


# vmin=6.
# vmax=17.
# levels=np.arange(6,16,0.5)

#======================================================================

ax1=plt.subplot(1, 1, 1)
# create Basemap instance for Robinson projection.
m = Basemap(projection='robin',lon_0=0,resolution='l') #globe


# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)
#m.shadedrelief()
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
#
#ax1 = m.contourf(x,y,DIFF_MAP,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
ax1 = m.contourf(x,y,CAM_MAP2,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)

#m.shadedrelief()
# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
#cbar.set_label(r'$1\sigma$, $\delta D_{WAX}$')
cbar.set_label(r'97.5% - 2.5% CI, $\delta D_{WAX}$')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'iCAM5 $\delta D_{C29}$, $\epsilon_{VEG}$ Monte Carlo 95% CI',fontsize=16)
#plt.title(r'iCAM5 $\delta D_{C29}$, $\epsilon_{VEG}$ Monte Carlo $1\sigma$',fontsize=16)
plt.show()
#======================================================================


