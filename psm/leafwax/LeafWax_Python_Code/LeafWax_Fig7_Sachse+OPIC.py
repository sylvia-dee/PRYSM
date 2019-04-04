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
# OIPC (Load Data)
#======================================================================

dir='/rdf/sd75/sylvia/LEAFWAX/'

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
OIPC_dD_ANN=OIPC[12,:,:]
lon=OIPC.getLongitude()
lat=OIPC.getLatitude()
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

#======================================================================
# EPSILON VALUES
#======================================================================


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
# calculate actual dD wax value

# 1. OIPC EPS

#OIPC_MAP1 = ((OIPC_r + 1000.) * (EPS/1000. + 1.)) - 1000.
OIPC_MAP1 = ((OIPC_r + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

#======================================================================

# SCATTER PLOT! 

# Create Location-lists of data for each site in the network.
cd /rdf/sd75/sylvia/LEAFWAX/

import csv
from numpy import genfromtxt

data = genfromtxt('sachse.csv',
				delimiter = ',',
				names=True,
				dtype=None
				#mask module may be used when better understood
				#,usemask=True
				)
dD_scatter=data['dD_c29_avg']
Slats=data['lat']
Slons=data['lon']
std=data['st_dev']

#======================================================================
# MAPS // PLOTTING
#======================================================================

# # READ IN LAND MASK
# dir = '/disk/staff/sylvia/LEAFWAX/'
# file_name=dir + 'landmask.f.e13.FiPIPDC5CN.f19_f19.ctrl.cam.h0.1850-2014.prec.nc'

# f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

# cdms2.axis.latitude_aliases.append("Y") 
# cdms2.axis.longitude_aliases.append("X")
# cdms2.axis.time_aliases.append("T")

# f.listdimension()
# f.listvariables()

# landmask = f('landmask',longitude=(0,360), latitude = (-90., 90.))

# thresh=0.5
# land_threshold_indices = landmask > thresh
# ocean_threshold_indices = landmask < thresh

# landmask[land_threshold_indices] = 1.0
# landmask[ocean_threshold_indices] = 0.0

# CAM_WAX_land_xy=CAM_WAX_xy*landmask
# CAM_WAX_land_pr=CAM_WAX_pr*landmask

# CAM_XY_mask=np.ma.masked_inside(CAM_WAX_land_xy, -0.1, 0.1)
# CAM_PR_mask=np.ma.masked_inside(CAM_WAX_land_pr, -0.1, 0.1)

#======================================================================
#======================================================================
from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

C_lat=latNp
C_lon=lonNp

MAP1,C_lon2 = shiftgrid(180.,OIPC_MAP1,C_lon,start=False)

# vmin=-400.
# vmax=-75.
# levels=np.arange(-450,-70,10)


vmin=-330.
vmax=0.
levels=np.arange(-380.,0,10.)


#======================================================================

ax1=plt.subplot(2, 1, 1)
# create Basemap instance for Robinson projection.
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=179, urcrnrlat=89)


# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
#m.drawcountries(linewidth=1.5)
#m.shadedrelief()
# compute the lons and lats to fit the projection
# x, y = m(*np.meshgrid(C_lon2, C_lat))

# # draw filled contours.
# # cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
# #
# #ax1 = m.contourf(x,y,DIFF_MAP,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
# ax1 = m.contourf(x,y,MAP1,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax,label='iCAM5-Xylem Water')

# plot scatter of data on top of map

vals = dD_scatter

a, b = m(Slons, Slats)
ax3=m.scatter(a,b,c=vals,cmap=plt.cm.nipy_spectral,marker='p',s=75,vmin=vmin,vmax=vmax,edgecolors='k',label='Observations')
plt.show()
#plt.legend()
#m.shadedrelief()
# add colorbar.
#cbar = m.colorbar(ax1,location='bottom',pad="2%")
#cbar.set_label(r'$\delta D_{WAX}$-PR')
#cbar.ax.tick_params(labelsize=16) 

plt.title(r'A. OIPC $\delta D_{C29}$, Precipitation',fontsize=16)

#======================================================================

ax1=plt.subplot(2, 1, 2)
# create Basemap instance for Robinson projection.
#m = Basemap(projection='robin',lon_0=0,resolution='l') #globe
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-65, urcrnrlon=180, urcrnrlat=89)


# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
#m.drawcountries(linewidth=1.5)
#m.shadedrelief()
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, C_lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
#
#ax1 = m.contourf(x,y,DIFF_MAP,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
ax1 = m.contourf(x,y,SPE_MAP1,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax,label='iCAM5-Xylem Water')

# plot scatter of data on top of map

vals = dD_scatter

a, b = m(Slons, Slats)
ax3=m.scatter(a,b,c=vals,cmap=plt.cm.nipy_spectral,marker='p',s=75,vmin=vmin,vmax=vmax,edgecolors='k',label='Observations')

#plt.legend()
#m.shadedrelief()
# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$-PR')
cbar.ax.tick_params(labelsize=16) 

plt.title(r'C. SPEEDY-IER $\delta D_{C29}$, Precipitation',fontsize=16)
plt.show()

#======================================================================

