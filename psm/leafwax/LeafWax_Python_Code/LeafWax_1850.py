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

start_time=cdtime.comptime(1980)
print 'start_time = ', start_time
end_time=cdtime.comptime(2005)
print 'end_time =', end_time

#======================================================================
# CAM5
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

#======================================================================

# # SAVE

# cd /home/geovault-02/sdee/LEAFWAX/

# np.save('CAM5_PR.npy',np.array(CAM5_PR))
# np.save('CAM5_dD_pr.npy',np.array(CAM5_dD_pr))
# np.save('CAM5_dD_xy.npy',np.array(CAM5_dD_xy))


#======================================================================
# SPEEDY
#======================================================================


# 1. Load input data/all environmental variables

dir = '/disk/staff/sylvia/SPEEDY/SPEEDY-IER_10/output/exp_R80/'
file_name=dir + 'attmR80.ctl'

#append file name

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

PRLH20 = f('prlh2o',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))
PRLHDO = f('prlhdo',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))
PRCH20 = f('prch2o',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))
PRCHDO = f('prchdo',longitude=(0,360), latitude = (-90., 90.), time=(start_time,end_time))


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

# Check dimensions

SPEEDY_dP.shape

# np.save('SPEEDY_dD_pr.npy',np.array(SPEEDY_dP))
# np.save('SPEEDY_dD_sw.npy',np.array(SPEEDY_dS))

#======================================================================
# 2.0 CALCULATE GROWING SEASON/ANNUAL/SEASONAL dDp, P
#======================================================================

# SPEEDY dD of Precip
SPEEDY_dD_pr_ANN=cdutil.YEAR(SPEEDY_dP)

# CAM5 dD of Precip

CAM_dD_pr_ANN=cdutil.YEAR(CAM5_dD_pr)

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

SPEEDY_pr_r=SPEEDY_dD_pr_ANN.regrid(grid3, regridTool='esmf', regridMethod='conserve')

#======================================================================
# EPSILON VALUES
#======================================================================

#======================================================================
# 2. PFT VEGETATION: READ IN FROM ETOT
#======================================================================

dir = '/disk/staff/sylvia/LEAFWAX/'

file_name=dir + 'Etot.year.2000.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#file_name=dir +'Etot.year.1850.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#append file name

#file_name=dir+'Etot.year.-19000.trace21k.pft.lgm.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

ETOT_2000 = f('Etot',longitude=(0,360), latitude = (-90., 90.))
# ETOT_1850= f('Etot',longitude=(0,360), latitude = (-90., 90.))

# create mask where Etot=0
ETOT_2000_mask=np.ma.masked_equal(ETOT_2000, -999.)
# SD seems to already be masked coming out of NCL.
#======================================================================

dir = '/disk/staff/sylvia/LEAFWAX/'

#file_name=dir + 'Etot.year.2000.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
file_name=dir +'Etot.year.1850.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#append file name

#file_name=dir+'Etot.year.-19000.trace21k.pft.lgm.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

#ETOT_2000 = f('Etot',longitude=(0,360), latitude = (-90., 90.))
ETOT_1850= f('Etot',longitude=(0,360), latitude = (-90., 90.))

# create mask where Etot=0
ETOT_1850_mask=np.ma.masked_equal(ETOT_1850, -999.)
# SD seems to already be masked coming out of NCL.
#======================================================================
# SET UP CALCULATION, PICK VARIABLE 
#======================================================================

# SET TO MAP AVERAGE

CAM_5_mask=np.ma.masked_greater(CAM5_dD_pr,1e+9)
CAM5_mask=CAM_5_mask.reshape(25,12,96,145)
CAM5_ANN_pr=np.mean(CAM5_mask,axis=(0,1))


# FIG 3/4
SPE_ANN=np.mean(SPEEDY_pr_r,axis=0)

#3. USE ETOT/EVEG (Figs 4,5 )

CAM_WAX_1850=((CAM5_ANN_pr + 1000.) * (ETOT_1850_mask/1000. + 1.)) - 1000.
SPE_WAX_1850=((SPE_ANN + 1000.) * (ETOT_1850_mask/1000. + 1.)) - 1000.

CAM_WAX_2000=((CAM5_ANN_pr + 1000.) * (ETOT_2000_mask/1000. + 1.)) - 1000.
SPE_WAX_2000=((SPE_ANN + 1000.) * (ETOT_2000_mask/1000. + 1.)) - 1000.

DIFFC=CAM_WAX_2000-CAM_WAX_1850
DIFFS=SPE_WAX_2000-SPE_WAX_1850

#======================================================================
#np.save('')
#======================================================================
# MAPS
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

# CAM_WAX_land=CAM_WAX*landmask
# SPE_WAX_land=SPE_WAX*landmask

DIFFC_land=DIFFC*landmask
DIFFS_land=DIFFS*landmask

CAM_MAP_mask=np.ma.masked_inside(DIFFC_land, -0.0001, 0.0001)
SPE_MAP_mask=np.ma.masked_inside(DIFFS_land, -0.0001, 0.0001)

#======================================================================

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica Neue']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

#======================================================================
from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

CAM_MAP2,C_lon2 = shiftgrid(180.,CAM_MAP_mask,C_lon,start=False)
SPE_MAP2,C_lon2 = shiftgrid(180.,SPE_MAP_mask,C_lon,start=False)

#MAPNOW=np.load('CAM5_1850_diff.npy')
vmin=-7.
vmax=7.
levels=np.arange(-14,14,0.5)

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

ax1 = m.contourf(x,y,CAM_MAP2,100,cmap=plt.cm.RdBu_r,levels=levels,vmin=vmin,vmax=vmax)

#m.shadedrelief()
# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ 2000-1850 Difference')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'iCAM5 Difference: 2000-1850 $\epsilon_{VEG}$',fontsize=16)
plt.show()
#======================================================================
