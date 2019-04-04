
#======================================================================

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
Sachse=data['dD_c29_avg']
Slats=data['lat']
Slons=data['lon']
std=data['st_dev']
#======================================================================

# Find points where both overlap -- nearest neighbor model point

ind_lat=np.zeros(Slats.shape)
for i in range(len(Slats)):
	latdiff=Slats[i]-C_lat
	abV=abs(latdiff)
	ind_lat[i]=np.argmin(abV)

ind_lon=np.zeros(Slons.shape)
for i in range(len(Slons)):
	londiff=Slons[i]-C_lon2
	abV=abs(londiff)
	ind_lon[i]=np.argmin(abV)

#======================================================================

Model_Data_xy=np.zeros(54)
Model_Data_pr=np.zeros(54)

for i in range(len(Slons)):
	Model_Data_xy[i]=CAM_MAP1[int(ind_lat[i]),int(ind_lon[i])]
	Model_Data_pr[i]=CAM_MAP2[int(ind_lat[i]),int(ind_lon[i])]

#======================================================================
# Compute RMSE
#======================================================================

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


# remove NaNs from XYLEM water dataset
PR=np.ma.masked_invalid(Model_Data_pr)
XY=np.ma.masked_invalid(Model_Data_xy)

RMS_pr=rmse(PR,Sachse)
RMS_xy=rmse(XY,Sachse)

# In [8]: RMS_pr
# Out[8]: 25.386984302781176

# In [9]: RMS_xy
# Out[9]: 29.891891689340802

# Create Map of Point Based Differences

DIFF_xy=XY-Sachse
DIFF_pr=PR-Sachse


#======================================================================
# USE REANALYSIS DATA TO CALCULATE E-P
#======================================================================

dir='/rdf/sd75/sylvia/LEAFWAX/'

file_name=dir + 'pevpr.mon.ltm.nc'

#append file name

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

# MANUALLY PULL DATES
EVAP = f('pevpr',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))
#W/m^2

#======================================================================

file_name=dir + 'prate.mon.ltm.nc'

f=cdms2.open(file_name) # NO QUOTES AROUND FILE NAME WHEN YOU OPEN IT!!!

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables()

PRECIP = f('prate',longitude=(0,360), latitude = (-90., 90.))#, time=(start_time,end_time))
#units: kg/m^2/s

NCEP_lons=np.array(PRECIP.getLongitude())
NCEP_lats=np.array(PRECIP.getLatitude())
#======================================================================

# convert W/m2 to mm/day and kg/m2/s to mm/day
# 1 Watt /m2 = 0.0864 MJ /m2/day
# 1 MJ /m2/day  =0.408 mm /day .

PRECIP_r=PRECIP*86400.
EVAP_r=EVAP*(0.0864)*(0.408)

PDIFF=np.mean(EVAP_r,axis=0) - np.mean(PRECIP_r,axis=0)

#======================================================================

# Find points where both overlap -- nearest neighbor model point

Slons2=Slons+180.

Nind_lat=np.zeros(Slats.shape)
for i in range(len(Slats)):
	latdiff=Slats[i]-NCEP_lats
	abV=abs(latdiff)
	Nind_lat[i]=np.argmin(abV)

Nind_lon=np.zeros(Slons2.shape)
for i in range(len(Slons)):
	londiff=Slons2[i]-NCEP_lons
	abV=abs(londiff)
	Nind_lon[i]=np.argmin(abV)

#======================================================================

Model_Data_PE=np.zeros(54)

for i in range(len(Slons)):
	Model_Data_PE[i]=PDIFF[int(Nind_lat[i]),int(Nind_lon[i])]

#======================================================================

# # numpy array -> axis -> grid
# lonNp = np.linspace(-180.0, 180.0, num=145)
# lonAxis = cdms2.createAxis(lonNp)
# latNp = np.linspace(-90., 90.0, num=96)
# latAxis = cdms2.createAxis(latNp)
# #======================================================================
# grid3 = cdms2.createRectGrid(latAxis, lonAxis, 'yx', type="generic")
# print grid3
# #======================================================================

# # HAVE TO USE THE ORIGINAL CDMS VARIABLE HERE, THEN RE-SAVE IF NECESSARY:

# PRECIP_r1 = PRECIP.regrid(grid3, regridTool='esmf', regridMethod='conserve')
# EVAP_r1 = EVAP.regrid(grid3, regridTool='esmf', regridMethod='conserve')
# lon=EVAP_r1.getLongitude()
# lat=EVAP_r1.getLatitude()

#======================================================================

plt.style.use('ggplot')

plt.subplot(1,1,1)
plt.scatter(Model_Data_PE,DIFF_pr,c='Navy',alpha=.5, s=200,marker='^',label='Precipitation-Forced')
plt.scatter(Model_Data_PE,DIFF_xy,c='orangered',alpha=.5, s=200,marker='o',label='Xylem-Forced')

# add colorbar.
plt.ylabel(r'$\delta D_{WAX}$ Residuals',fontsize=18)
plt.xlabel(r'E-P [mm/day] (NCEP REANALYSIS)',fontsize=18)
plt.tick_params(labelsize=16) 
plt.title(r'Model-Data Comparison, Impact of E-P',fontsize=20)
plt.ylim([-45,45])
plt.legend()


plt.show()
#======================================================================

