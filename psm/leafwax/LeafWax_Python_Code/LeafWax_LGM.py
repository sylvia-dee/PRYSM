
dir = '/disk/staff/sylvia/LEAFWAX/'

#file_name=dir + 'Etot.year.2000.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#file_name=dir +'Etot.year.1850.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#append file name

file_name=dir+'Etot.year.-19000.trace21k.pft.lgm.nc'

# LGM: the land mask changes between the two time periods, so a difference map of predicted 
# dDc29 LGM-2000 will probably give you errors. For this figure, I'd recommend just masking out 
# anything that has a missing value either during the LGM or during the year 2000. That will take c
# are of points like Sundaland that have data during the LGM but not 2000, and points like under 
# the Laurentide ice sheet that have data during 2000 but not during the LGM. 

# TWO VARIABLES:

# CAM_dD_pr_ANN
# SPEEDY_dD_pr_ANN


f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

ETOT = f('Etot',longitude=(0.,360.), latitude = (-90., 90.))
# ETOT_1850= f('Etot',longitude=(0,360), latitude = (-90., 90.))

# create mask where Etot=0
ETOT_mask=np.ma.masked_equal(ETOT, -999.)
# SD seems to already be masked coming out of NCL.


#======================================================================
# numpy array -> axis -> grid
lonNp = np.linspace(0.0, 360.0, num=97.)
lonAxis = cdms2.createAxis(lonNp)
latNp = np.linspace(-90., 90.0, num=48.)
latAxis = cdms2.createAxis(latNp)
#======================================================================
grid3 = cdms2.createRectGrid(latAxis, lonAxis, 'yx', type="generic")
print grid3
#======================================================================

# HAVE TO USE THE ORIGINAL CDMS VARIABLE HERE, THEN RE-SAVE IF NECESSARY:

CAM_dD_pr_regrid=CAM_dD_pr_ANN.regrid(grid3, regridTool='esmf', regridMethod='conserve')

#======================================================================
# SET UP CALCULATION, PICK VARIABLE 
#======================================================================

# SET TO MAP AVERAGE

CAM_5_mask=np.ma.masked_greater(CAM_dD_pr_regrid,1e+9)
CAM5_ANN_pr=np.mean(CAM_5_mask,axis=(0))

#======================================================================

#3. USE ETOT/EVEG (Figs 4,5 )

CAM_WAX_LGM=((CAM5_ANN_pr + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.
# SPE_WAX_LGM=((SPE_ANN + 1000.) * (ETOT_mask/1000. + 1.)) - 1000.

#======================================================================
# FIG 5 USE 2000 VALUES
#======================================================================

dir = '/disk/staff/sylvia/LEAFWAX/'

file_name=dir + 'Etot.year.2000.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'
#file_name=dir +'Etot.year.1850.surfdata.pftdyn_1.9x2.5_rcp4.5_simyr1850-2100_c100322.nc'

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

ETOT_2000 = f('Etot',longitude=(0.,360.), latitude = (-90., 90.))
# ETOT_1850= f('Etot',longitude=(0,360), latitude = (-90., 90.))

#======================================================================

lonNp = np.linspace(0.0, 360.0, num=97.)
lonAxis = cdms2.createAxis(lonNp)
latNp = np.linspace(-90., 90.0, num=48.)
latAxis = cdms2.createAxis(latNp)
#======================================================================
grid3 = cdms2.createRectGrid(latAxis, lonAxis, 'yx', type="generic")
print grid3
#======================================================================

# HAVE TO USE THE ORIGINAL CDMS VARIABLE HERE, THEN RE-SAVE IF NECESSARY:

ETOT_2000_regrid=ETOT_2000.regrid(grid3, regridTool='esmf', regridMethod='conserve')
#======================================================================
# create mask where Etot=0
ETOT_mask_2000=np.ma.masked_equal(ETOT_2000_regrid, -999.)
#======================================================================

CAM_WAX_2000=((CAM5_ANN_pr + 1000.) * (ETOT_mask_2000/1000. + 1.)) - 1000.
#SPE_WAX_2000=((SPE_ANN + 1000.) * (ETOT_mask_2000/1000. + 1.)) - 1000.

#======================================================================


DIFF_CAM=CAM_WAX_2000-CAM_WAX_LGM
#DIFF_SPE=SPE_WAX_2000-SPE_WAX_LGM

#======================================================================
# MAPS
#======================================================================

# READ IN LAND MASK
dir = '/disk/staff/sylvia/LEAFWAX/'
file_name=dir + 'landmask.f.e13.FiPIPDC5CN.f19_f19.ctrl.cam.h0.1850-2014.prec.nc'
f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

landmask_2000 = f('landmask',longitude=(0.,360.), latitude = (-90., 90.))
landmask_2000_regrid=landmask_2000.regrid(grid3, regridTool='esmf', regridMethod='conserve')

#======================================================================

file_name=dir+'landmask_LGM_TraCE.nc'
f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

landmask_LGM = f('landmask',longitude=(0.,360.), latitude = (-90., 90.))

#======================================================================

# thresh=0.5
# land_threshold_indices = landmask > thresh
# ocean_threshold_indices = landmask < thresh

# landmask[land_threshold_indices] = 1.0
# landmask[ocean_threshold_indices] = 0.0

# CAM_WAX_land=CAM_WAX_LGM*landmask
# SPE_WAX_land=SPE_WAX_LGM*landmask

# CAM_MAP_mask=np.ma.masked_inside(CAM_WAX_land, -0.1, 0.1)
# SPE_MAP_mask=np.ma.masked_inside(SPE_WAX_land, -0.1, 0.1)


# USE THIS FOR THE DIFFERENCE MAPS EMPLOYING BOTH LAND MASKS:

landmask_full=landmask_LGM+landmask_2000_regrid

thresh=1.5
land_threshold_indices = landmask_full > thresh
ocean_threshold_indices = landmask_full < thresh

landmask_full[land_threshold_indices] = 1.0
landmask_full[ocean_threshold_indices] = 0.0

#======================================================================
# Difference Maps, 2000-LGM

DIFF_CAM_WAX_land=DIFF_CAM*landmask_full
DIFF_CAM_MAP_mask=np.ma.masked_inside(DIFF_CAM_WAX_land, -0.2, 0.2)

#======================================================================
from mpl_toolkits.basemap import Basemap, shiftgrid
plt.style.use('ggplot')

lon=landmask_LGM.getLongitude()
lat=landmask_LGM.getLatitude()

lon=np.array(lon)
lat=np.array(lat)

lon=lon[0:96]

CAM_MAP2,C_lon2 = shiftgrid(180.,DIFF_CAM_MAP_mask[:,0:96],lon,start=False)

vmin=-60.
vmax=60.
levels=np.arange(-90,90,1)

#======================================================================

ax1=plt.subplot(1, 1, 1)
# create Basemap instance for Robinson projection.
m = Basemap(projection='robin',lon_0=0,resolution='l') #globe


# draw costlines and coutries
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.5)
#m.shadedrelief()
# compute the lons and lats to fit the projection
x, y = m(*np.meshgrid(C_lon2, lat))

# draw filled contours.
# cs = m.contour(x,y,d18o,13,colors='black',linewidths=0.25)
#
#ax1 = m.contourf(x,y,DIFF_MAP,100,cmap=plt.cm.nipy_spectral,levels=levels,vmin=vmin,vmax=vmax)
ax1 = m.contourf(x,y,CAM_MAP2,100,cmap=plt.cm.RdBu_r,levels=levels,vmin=vmin,vmax=vmax)

#m.shadedrelief()
# add colorbar.
cbar = m.colorbar(ax1,location='bottom',pad="2%")
cbar.set_label(r'$\delta D_{WAX}$ (2000-LGM $\epsilon_{VEG}$)')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'iCAM5: 2000 - LGM $\epsilon_{VEG}$',fontsize=16)
plt.show()
#======================================================================



ax2=plt.subplot(1, 2, 1)
m = Basemap(projection='robin',lon_0=0,resolution='l') #globe

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
plt.title(r'OIPC DJF $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)



#======================================================================
#======================================================================

ax2=plt.subplot(1, 2, 2)
m = Basemap(projection='robin',lon_0=0,resolution='l') #globe

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
cbar.set_label(r'$\delta D_P$ JJA')
cbar.ax.tick_params(labelsize=16) 
plt.title(r'OIPC JJA $\delta D_{WAX}$, Fixed $\epsilon = -120$',fontsize=16)

plt.show()

#======================================================================