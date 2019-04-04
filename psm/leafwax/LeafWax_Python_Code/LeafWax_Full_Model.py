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
import sys, os, cdtime, cdutil, cdms2,  MV2, time, datetime
import numpy as np
import matplotlib.pyplot as plt
from math import exp

#======================================================================
# 1.0 Load Data for each Climate Model
#======================================================================

start_time=cdtime.comptime(1980)
print 'start_time = ', start_time
end_time=cdtime.comptime(2005)
print 'end_time =', end_time


#======================================================================
# OIPC
#======================================================================

dir='/rdf/sd75/sylvia/LEAFWAX/'

file_name=dir + 'oipc_global10_v2.nc'

#append file name

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

# MANUALLY PULL DATES

OIPC = f('oipc_global10_v2',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))


#======================================================================
# CAM5
#======================================================================

dir='/rdf/sd75/sylvia/LEAFWAX/'

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

# # SAVE

# cd /home/geovault-02/sdee/LEAFWAX/

# np.save('CAM5_PR.npy',np.array(CAM5_PR))
# np.save('CAM5_dD_pr.npy',np.array(CAM5_dD_pr))
# np.save('CAM5_dD_xy.npy',np.array(CAM5_dD_xy))


#======================================================================
# SPEEDY
#======================================================================

# PRELOAD from directory:
# SPEEDY_dD_pr_ANN=np.load('SPEEDY_dD_pr.npy')

# # 1. Load input data/all environmental variables

dir='/disk/staff/sylvia/SPEEDY/SPEEDY-IER_10/output/exp_R80/'
file_name=dir + 'attmR80.ctl'

#append file name

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

PRLH20 = f('prlh2o',longitude=(0,360), latitude = (-90., 90.),time=(start_time,end_time))
PRLHDO = f('prlhdo',longitude=(0,360), latitude = (-90., 90.),time=(start_time,end_time))
PRCH20 = f('prch2o',longitude=(0,360), latitude = (-90., 90.),time=(start_time,end_time))
PRCHDO = f('prchdo',longitude=(0,360), latitude = (-90., 90.),time=(start_time,end_time))


lat = PRLH20.getLatitude()
lon = PRLH20.getLongitude()

cdutil.setTimeBoundsMonthly(PRLH20)
cdutil.setTimeBoundsMonthly(PRLHDO)
cdutil.setTimeBoundsMonthly(PRCH20)
cdutil.setTimeBoundsMonthly(PRCHDO)

IH20 = PRLH20 + PRCH20
IHDO = PRLHDO + PRCHDO


# Check dimension

IH20.shape
IHDO.shape

SPEEDY_dP=((IHDO/IH20)-1)*1000

#  Soil Water

# WTWBH2O = f('wtwbh2o',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))
# WTWBDO = f('wtwbdo',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))

# cdutil.setTimeBoundsMonthly(WTWBH2O)
# cdutil.setTimeBoundsMonthly(WTWBDO)

# Soil_R18=(WTWBDO/WTWBH2O)
# SPEEDY_dS = (Soil_R18-1)*1000

# Check dimensions

SPEEDY_dP.shape
# SPEEDY_dS.shape

# np.save('SPEEDY_dD_pr.npy',np.array(SPEEDY_dP))
# np.save('SPEEDY_dD_sw.npy',np.array(SPEEDY_dS))

#======================================================================
# 2.0 CALCULATE GROWING SEASON/ANNUAL/SEASONAL dDp, P
#======================================================================

# SPEEDY dD of Precip
SPEEDY_dD_pr_ANN=cdutil.YEAR(SPEEDY_dP)
# SPEEDY_dD_pr_DJF=cdutil.DJF(SPEEDY_dP)
# SPEEDY_dD_pr_MAM=cdutil.MAM(SPEEDY_dP)
# SPEEDY_dD_pr_JJA=cdutil.JJA(SPEEDY_dP)
# SPEEDY_dD_pr_SON=cdutil.SON(SPEEDY_dP)

# # SPEEDY dD of Soil Moisture
# SPEEDY_dD_sw_ANN=cdutil.YEAR(SPEEDY_dS)
# SPEEDY_dD_sw_DJF=cdutil.DJF(SPEEDY_dS)
# SPEEDY_dD_sw_MAM=cdutil.MAM(SPEEDY_dS)
# SPEEDY_dD_sw_JJA=cdutil.JJA(SPEEDY_dS)
# SPEEDY_dD_sw_SON=cdutil.SON(SPEEDY_dS)

OIPC_dD_ANN=OIPC[12,:,:]
OIPC_dD_DJF=OIPC[0,:,:]
OIPC_dD_MAM=OIPC[2:4,:,:]
OIPC_dD_JJA=OIPC[6,:,:]
OIPC_dD_SON=OIPC[8:10,:,:]

# CAM5 dD of Precip

CAM_dD_pr_ANN=cdutil.YEAR(CAM5_dD_pr)
# CAM_dD_pr_DJF=cdutil.DJF(CAM5_dD_pr)
# CAM_dD_pr_MAM=cdutil.MAM(CAM5_dD_pr)
# CAM_dD_pr_JJA=cdutil.JJA(CAM5_dD_pr)
# CAM_dD_pr_SON=cdutil.SON(CAM5_dD_pr)

# CAM5 dD of Xylem Water
CAM_dD_xy_ANN=cdutil.YEAR(CAM5_dD_xy)
# CAM_dD_xy_DJF=cdutil.DJF(CAM5_dD_xy)
# CAM_dD_xy_MAM=cdutil.MAM(CAM5_dD_xy)
# CAM_dD_xy_JJA=cdutil.JJA(CAM5_dD_xy)
# CAM_dD_xy_SON=cdutil.SON(CAM5_dD_xy)

# CAM5 Precip Amount
CAM_pr_ANN=cdutil.YEAR(CAM5_PR)
# CAM_pr_DJF=cdutil.DJF(CAM5_PR)
# CAM_pr_MAM=cdutil.MAM(CAM5_PR)
# CAM_pr_JJA=cdutil.JJA(CAM5_PR)
# CAM_pr_SON=cdutil.SON(CAM5_PR)

#======================================================================
# NEEDED?
#; Transform P from m/s to mm/day
#P = P * 86400000. 
#======================================================================
#======================================================================
# 3. REGRID ALL DATA TO CAM GRID
#======================================================================
#======================================================================
#======================================================================
# numpy array -> axis -> grid
lonNp = np.linspace(0.0, 360.0, num=145)
lonAxis = cdms2.createAxis(lonNp)
latNp = np.linspace(-90., 90.0, num=96)
latAxis = cdms2.createAxis(latNp)
#======================================================================
grid3 = cdms2.createRectGrid(latAxis, lonAxis, 'yx', type="generic")
print grid3
#======================================================================

# HAVE TO USE THE ORIGINAL CDMS VARIABLE HERE, THEN RE-SAVE IF NECESSARY:

OIPC_r = OIPC_dD_ANN.regrid(grid3, regridTool='esmf', regridMethod='conserve')
OIPC_DJF_r=OIPC_dD_DJF.regrid(grid3, regridTool='esmf', regridMethod='conserve')
OIPC_JJA_r=OIPC_dD_JJA.regrid(grid3, regridTool='esmf', regridMethod='conserve')


SPEEDY_pr_r=SPEEDY_dD_pr_ANN.regrid(grid3, regridTool='esmf', regridMethod='conserve')
# SPEEDY_sw_r=SPEEDY_dD_sw_ANN.regrid(grid3, regridTool='esmf', regridMethod='conserve')

# SPEEDY_pr_DJF_r=SPEEDY_dD_pr_DJF.regrid(grid3, regridTool='esmf', regridMethod='conserve')
# SPEEDY_pr_JJA_r=SPEEDY_dD_pr_JJA.regrid(grid3, regridTool='esmf', regridMethod='conserve')

#======================================================================
# EPSILON VALUES
#======================================================================

# 1. DONT USE ETOT: just use an average epsilon value.

EPS=-120.


#======================================================================
# 2. PFT VEGETATION: READ IN FROM ETOT
#======================================================================

dir='/rdf/sd75/sylvia/LEAFWAX/'

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
# SET UP CALCULATION, PICK VARIABLE 
#======================================================================

# SET TO MAP AVERAGE

CAM_5_mask=np.ma.masked_greater(CAM5_dD_pr,1e+9)
CAM5_mask=CAM_5_mask.reshape(25,12,96,145)
CAM5_ANN_pr=np.mean(CAM5_mask,axis=(0,1))

CAM_5_mask=np.ma.masked_greater(CAM5_dD_xy,1e+9)
CAM5_mask=CAM_5_mask.reshape(25,12,96,145)
CAM5_ANN_xy=np.mean(CAM5_mask,axis=(0,1))


# FIG 6
# CAM5_ANN_xy=np.mean(CAM_5_mask,axis=0)
# SPE_ANN_sw=np.mean(SPEEDY_sw_r,axis=0)
#CAM_5_mask=np.ma.masked_greater(CAM_dD_xy_ANN,1e+9)
#CAM5_mask=CAM_5_mask.reshape(25,12,96,145)


# FIG 3/4
SPE_ANN=np.mean(SPEEDY_pr_r,axis=0)

# VAR1=np.average(CAM_dD_xy_ANN,axis=0)
# #VAR1_mask=np.ma.masked_equal(VAR1, 0.)

# #VAR2=np.average(CAM_dD_pr_ANN,axis=0)
# VAR2=np.ma.average(CAM5_ANN,axis=0)


#======================================================================
# for sesitivity tests: dD precip, dD soil water, 
# dD xylem water, ANN, DJF, MAM, JJA, SON
#======================================================================

# FIXED EPS

# 1. OIPC EPS

OIPC_MAP1 = ((OIPC_r + 1000.) * (EPS/1000. + 1.)) - 1000.
OIPC_MAPDJF=((OIPC_DJF_r + 1000.) * (EPS/1000. + 1.)) - 1000.
OIPC_MAPJJA=((OIPC_JJA_r + 1000.) * (EPS/1000. + 1.)) - 1000.

#2. CAM, SPEEDY EPS

CAM_WAX=((CAM5_ANN_pr + 1000.) * (EPS/1000. + 1.)) - 1000.
CAM_WAX_XY=((CAM5_ANN_xy + 1000.) * (EPS/1000. + 1.)) - 1000.


SPE_WAX=((SPE_ANN + 1000.) * (EPS/1000. + 1.)) - 1000.

# SAVE FOR SCATTERPLOTS

np.save('OIPC_ANN_FIXED_EPS.npy',np.array(OIPC_MAP1))
np.save('CAM_PR_ANN_FIXED_EPS.npy',np.array(CAM_WAX))
np.save('CAM_XY_ANN_FIXED_EPS.npy',np.array(CAM_WAX_XY))
np.save('SPEEDY_ANN_FIXED_EPS.npy',np.array(SPE_WAX))

#======================================================================
#3. USE ETOT/EVEG (based on plant functional type)
#======================================================================

OIPC_MAP1 = ((OIPC_r + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

CAM_WAX=((CAM5_ANN + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

SPE_WAX=((SPE_ANN + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

#======================================================================
# FIG 6 USE dD XY and SW
#======================================================================

# CAM_WAX_xy=((CAM5_ANN_xy + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.
# SPE_WAX_xy=((SPE_ANN_sw + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

# FIG 7:
CAM_WAX_pr=((CAM5_ANN_pr + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.
CAM_WAX_xy=((CAM5_ANN_xy + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

DIFF=CAM_WAX_xy-CAM_WAX_pr

# dDwax_PR = ((VAR2 + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

# SAVE FOR SCATTERPLOTS

np.save('OIPC_ANN_ETOT_EPS.npy',np.array(OIPC_MAP1))
np.save('CAM_PR_ANN_ETOT_EPS.npy',np.array(CAM_WAX_pr))
np.save('CAM_XY_ANN_ETOT_EPS.npy',np.array(CAM_WAX_xy))
np.save('SPEEDY_ANN_ETOT_EPS.npy',np.array(SPE_WAX))

#======================================================================
#======================================================================
# MAPS
#======================================================================
#======================================================================

C_lat=np.array(CAM5_PR.getLatitude())
C_lon=np.array(CAM5_PR.getLongitude())


# READ IN LAND MASK
dir = '/disk/staff/sylvia/LEAFWAX/'
file_name=dir + 'landmask.f.e13.FiPIPDC5CN.f19_f19.ctrl.cam.h0.1850-2014.prec.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

landmask = f('landmask',longitude=(0,360), latitude = (-90., 90.))

thresh=0.5
land_threshold_indices = landmask > thresh
ocean_threshold_indices = landmask < thresh

landmask[land_threshold_indices] = 1.0
landmask[ocean_threshold_indices] = 0.0

#======================================================================

OANN=OIPC_MAP1*landmask
ODJF=OIPC_MAPDJF*landmask
OJJA=OIPC_MAPJJA*landmask

#======================================================================

OIPC_MAP_mask2=np.ma.masked_inside(OANN, -0.1, 0.1)
OIPC_MAP_DJF_mask2=np.ma.masked_inside(ODJF, -0.1, 0.1)
OIPC_MAP_JJA_mask2=np.ma.masked_inside(OJJA, -0.1, 0.1)

#======================================================================

CAM_WAX_land=CAM_WAX*landmask
SPE_WAX_land=SPE_WAX*landmask

# Fig 7
# DIFF_land=DIFF*landmask

# CAM_WAX_land=CAM_WAX_xy*landmask
# SPE_WAX_land=SPE_WAX_xy*landmask

CAM_MAP_mask=np.ma.masked_inside(CAM_WAX_land, -0.1, 0.1)
SPE_MAP_mask=np.ma.masked_inside(SPE_WAX_land, -0.1, 0.1)

#DIFF_mask=np.ma.masked_inside(DIFF_land, -0.1, 0.1)
#======================================================================

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

#======================================================================
from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

#C1=np.average(CAM_dD_pr_ANN,axis=0)
#S1=np.average(SPEEDY_dD_pr_ANN,axis=0)

OIPC_MAP1,C_lon2 = shiftgrid(180.,OIPC_MAP_mask2,C_lon,start=False)
OIPC_MAP2,C_lon2 = shiftgrid(180.,OIPC_MAP_DJF_mask2,C_lon,start=False)
OIPC_MAP3,C_lon2 = shiftgrid(180.,OIPC_MAP_JJA_mask2,C_lon,start=False)

CAM_MAP2,C_lon2 = shiftgrid(180.,CAM_MAP_mask,C_lon,start=False)
SPE_MAP2,C_lon2 = shiftgrid(180.,SPE_MAP_mask,C_lon,start=False)

#DIFF_MAP,C_lon2 = shiftgrid(180.,DIFF_mask,C_lon,start=False)

# S2,S_lon2 = shiftgrid(180.,S1,S_lon,start=False)
# DWAX, C_lon2= shiftgrid(180.,dDwax_XY,C_lon,start=False)
# DWAX2, C_lon2= shiftgrid(180.,dDwax_PR,C_lon,start=False)

# Cr_value,newlons = shiftgrid(180.,Cr_value,xlons,start=False)

vmin=-330.
vmax=0.
levels=np.arange(-380.,0,10.)

#======================================================================
# MANUSCRIPT FIGURE 2
#======================================================================

plt.style.use('ggplot')

ax1=plt.subplot(1, 1, 1)
# create Basemap instance for Robinson projection.
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)


# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)
#m.shadedrelief()
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
#
ax1 = m.contourf(x,y,OIPC_MAP1,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
#ax1 = m.contourf(x,y,DIFF_MAP,100,cmap=plt.cm.PuOr_r,levels=levels,vmin=vmin,vmax=vmax)

#m.shadedrelief()
# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$, Annual Average')
cbar.ax.tick_params(labelsize=16) 
plt.title('A. OIPC (ANN) $\delta D_{WAX}$: Fixed $\epsilon = -120$',fontsize=16)
plt.show()
#======================================================================




ax2=plt.subplot(1, 2, 1)
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)

# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)

# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
ax1 = m.contourf(x,y,OIPC_MAP2,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
#m.shadedrelief()

# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ DJF')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'B. OIPC DJF $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)



#======================================================================
#======================================================================

ax2=plt.subplot(1, 2, 2)
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)

# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)

# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
ax1 = m.contourf(x,y,OIPC_MAP3,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)

# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ JJA')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'C. OIPC JJA $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)

plt.show()

#======================================================================

# ax2=plt.subplot(2, 2, 3)
# m = Basemap(projection='robin',lon_0=0,resolution='l') #globe

# # draw costlines and coutries
# m.drawcoastlines(linewidth=1.5)
# m.drawcountries(linewidth=1.5)

# # compute the lons and lats to fit the projection
# x, y = m(*np.meshgrid(I_lon2, I_lat))

# # draw filled contours.
# # cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
# ax1 = m.contourf(x,y,I2,50,cmap=plt.cm.RdBu_r,levels=levels)

# # add colorbar.
# cbar = m.colorbar(ax1,location='bottom',pad="2%")
# cbar.set_label(r'$\delta D_P$ ANNUAL AVERAGE')
# cbar.ax.tick_params(labelsize=16) 
# plt.title(r'IsoGSM (nudged) $\delta D_P$, Modern',fontsize=16)

# #======================================================================

# ax4=plt.subplot(2, 2, 4)
# m = Basemap(projection='robin',lon_0=0,resolution='l') #globe

# # draw costlines and coutries
# m.drawcoastlines(linewidth=1.5)
# m.drawcountries(linewidth=1.5)

# # compute the lons and lats to fit the projection
# x, y = m(*np.meshgrid(S_lon2, S_lat))

# # draw filled contours.
# # cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
# ax1 = m.contourf(x,y,S2,50,cmap=plt.cm.RdBu_r,levels=levels)

# # add colorbar.
# cbar = m.colorbar(ax1,location='bottom',pad="2%")
# cbar.set_label(r'$\delta D_P$ ANNUAL AVERAGE')
# cbar.ax.tick_params(labelsize=16) 
# plt.title(r'SPEEDY $\delta D_P$, Modern',fontsize=16)

# #======================================================================
# plt.suptitle(r'$\delta D$ of Precipitation: Multi-Model Comparison')
# plt.show()

#======================================================================
# MANUSCRIPT FIGURE 3
#======================================================================


ax2=plt.subplot(2, 1, 1)
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)

# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)

# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
ax1 = m.contourf(x,y,CAM_MAP2,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
#m.shadedrelief()

# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ (ANN)')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'A. iCAM5 $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)



#======================================================================
#======================================================================

ax2=plt.subplot(2, 1, 2)
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)

# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)

# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
ax1 = m.contourf(x,y,SPE_MAP2,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)

# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ (ANN)')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'B. SPEEDY $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)

plt.show()



#======================================================================
# MANUSCRIPT FIGURE 6
#======================================================================

CAM_WAX_pr=((CAM5_ANN_pr + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.
CAM_WAX_xy=((CAM5_ANN_xy + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

DIFF=CAM_WAX_xy-CAM_WAX_pr
#======================================================================

from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

DIFFMAP,C_lon2 = shiftgrid(180.,DIFF,C_lon,start=False)

vmin=-80.
vmax=80.
levels=np.arange(-300.,310,10.)
#======================================================================


ax2=plt.subplot(1, 1, 1)
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)

# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
#m.drawcountries(linewidth=1.5)

# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.

ax1 = m.contourf(x,y,DIFFMAP,100,cmap=plt.cm.BrBG_r,levels=levels,vmin=vmin,vmax=vmax)

# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ (XYLEM minus PRECIP)')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'iCAM5 Difference: $\delta D_{WAX-XYLEM}$ minus $\delta D_{WAX-PRECIP}$, $\epsilon_{VEG}$',fontsize=16)

plt.show()

