from pylab import *
import sys, os, cdtime, cdutil, cdms2,  MV2, time, datetime
import numpy as np
import matplotlib.pyplot as plt
from math import exp

#import Scientific.IO.NetCDF 

#======================================================================

# Load Datasets.
e#=1.1*np.power(5.,-0.73)
#======================================================================

# Use full grid for SODA
dir = '/Users/sdee/Desktop/Sediment Cores/ERA-Interim Tuning/'
file_name=dir + '_grib2netcdf-atls20-a562cefde8a29a7288fa0b8b7f9413f7-BXvP0w.nc'
f=cdms2.open(file_name)

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables() # ['t2m', 'v10', 'sp', 'd2m', 'si10', 'u10']

#2m Air#======================================================================

t2m_M= f('t2m',longitude=(33, 36), latitude = (-14, -9))
t2m_T= f('t2m',longitude=(28, 31), latitude = (-8, -3))

t2m_M.shape
t2m_T.shape

T_T2M = cdutil.averager(t2m_T, axis='xy')
M_T2M = cdutil.averager(t2m_M, axis='xy')

#2m Dew#======================================================================

d2m_M= f('d2m',longitude=(33, 36), latitude = (-14, -9))
d2m_T= f('d2m',longitude=(28, 31), latitude = (-8, -3))

d2m_M.shape
d2m_T.shape

T_D2M = cdutil.averager(d2m_T, axis='xy')
M_D2M = cdutil.averager(d2m_M, axis='xy')

# U-wind#======================================================================

u10_M= f('u10',longitude=(33, 36), latitude = (-14, -9))
u10_T= f('u10',longitude=(28, 31), latitude = (-8, -3))

u10_M.shape
u10_T.shape

T_u10 = cdutil.averager(u10_T, axis='xy')
M_u10 = cdutil.averager(u10_M, axis='xy')

# V-wind#======================================================================

v10_M= f('v10',longitude=(33, 36), latitude = (-14, -9))
v10_T= f('v10',longitude=(28, 31), latitude = (-8, -3))

v10_M.shape
v10_T.shape

T_v10 = cdutil.averager(v10_T, axis='xy')
M_v10 = cdutil.averager(v10_M, axis='xy')

# WindSpeed#======================================================================


T_WIND=np.sqrt(np.power(T_u10,2)+np.power(T_v10,2))
M_WIND=np.sqrt(np.power(M_u10,2)+np.power(M_v10,2))

#Surface Pressure#======================================================================

sp_M= f('sp',longitude=(33, 36), latitude = (-14, -9))
sp_T= f('sp',longitude=(28, 31), latitude = (-8, -3))

sp_M.shape
sp_T.shape

T_sp = cdutil.averager(sp_T, axis='xy')
M_sp = cdutil.averager(sp_M, axis='xy')

#Wind Speed 10m#======================================================================

si10_M= f('si10',longitude=(33, 36), latitude = (-14, -9))
si10_T= f('si10',longitude=(28, 31), latitude = (-8, -3))

si10_M.shape
si10_T.shape

T_si10 = cdutil.averager(si10_T, axis='xy')
M_si10 = cdutil.averager(si10_M, axis='xy')


#Radiation#======================================================================

dir = '/Users/sdee/Desktop/Sediment Cores/ERA-Interim Tuning/'
file_name=dir + '_grib2netcdf-atls00-a562cefde8a29a7288fa0b8b7f9413f7-D4ZcNr.nc'
f=cdms2.open(file_name)

cdms2.axis.latitude_aliases.append("Y") 
cdms2.axis.longitude_aliases.append("X")
cdms2.axis.time_aliases.append("T")

f.listdimension()
f.listvariables() # ['t2m', 'v10', 'sp', 'd2m', 'si10', 'u10']

ssrd_M= f('ssrd',longitude=(33, 36), latitude = (-14, -9))
ssrd_T= f('ssrd',longitude=(28, 31), latitude = (-8, -3))

ssrd_M.shape
ssrd_T.shape

T_ssrd = cdutil.averager(ssrd_T, axis='xy')
M_ssrd = cdutil.averager(ssrd_M, axis='xy')

#Radiation#======================================================================

strd_M= f('strd',longitude=(33, 36), latitude = (-14, -9))
strd_T= f('strd',longitude=(28, 31), latitude = (-8, -3))

strd_M.shape
strd_T.shape

T_strd = cdutil.averager(strd_T, axis='xy')
M_strd = cdutil.averager(strd_M, axis='xy')

#CONVERSIONS======================================================================
# Convert all variables to /month
# 1/s * 1min/60s * 1hr/60min * 1day/24 hours * 1 month*30 days

T_strd_mon=T_strd*(1./60.*1./60.*1./24.)#*1./30.)
M_strd_mon=M_strd*(1./60.*1./60.*1./24.)#*1./30.)
T_ssrd_mon=T_ssrd*(1./60.*1./60.*1./24.)#*1./30.)
M_ssrd_mon=M_ssrd*(1./60.*1./60.*1./24.)#*1./30.)

# Convert Dewpoint temperature to RH

T_RH = 100.-5.*(T_T2M-T_D2M)
M_RH = 100.-5.*(M_T2M-M_D2M)

# BIAS======================================================================
# BIAS CORRECT WIND FIELD AND RH (OPTIONAL)

T_RH.shape

RH_bias=np.array((-20.56,-19.05,-20.57,-17.66,-4.59,9.09, 17.95,25.22,22.4,31.53,-21.62,-26.81))
WIND_bias=np.array((3.36,2.41,2.18,2.89,2.58,2.62,1.42,0.6,1.19,2.2,2.79,3.24))

RHcorr=np.tile(RH_bias,37)
Wcorr=np.tile(WIND_bias,37)

T_RH_corr=T_RH+RHcorr
T_WIND_corr=T_WIND+Wcorr
#Time Axis#======================================================================
years=np.arange(1,38,1)
year_array=np.zeros((37,12))
for i in range(37):
	year_array[i,:]=years[i]


year_array=year_array.reshape(444)


days=np.zeros(459)
days[0]=15.
for i in range(len(days)-1):
	days[i+1]=days[i]+np.floor(365./12.)

# time=strd_M.getTime()
# time=np.array(time)

# time_month=1900+(time/24./365.)

#Save#======================================================================

ERA_INTERIM_1979_2016=np.zeros((444,8))
ERA_INTERIM_1979_2016[:,0]=year_array
ERA_INTERIM_1979_2016[:,1]=days[0:-15]
ERA_INTERIM_1979_2016[:,2]=np.array(T_T2M)
ERA_INTERIM_1979_2016[:,3]=np.array(T_RH_corr)
ERA_INTERIM_1979_2016[:,4]=np.array(T_WIND_corr)
ERA_INTERIM_1979_2016[:,5]=np.array(T_ssrd_mon[0:-15])
ERA_INTERIM_1979_2016[:,6]=np.array(T_strd_mon[0:-15])
ERA_INTERIM_1979_2016[:,7]=np.array(T_sp)

np.savetxt('ERA_INTERIM_1979_2016_Tanganyika_BIASCORRECT.txt',ERA_INTERIM_1979_2016,'%5.2f',delimiter='\t')


#Save Malawi#======================================================================

ERA_INTERIM_1979_2016=np.zeros((444,8))
ERA_INTERIM_1979_2016[:,0]=year_array
ERA_INTERIM_1979_2016[:,1]=days[0:-15]
ERA_INTERIM_1979_2016[:,2]=np.array(M_T2M)
ERA_INTERIM_1979_2016[:,3]=np.array(M_RH)
ERA_INTERIM_1979_2016[:,4]=np.array(M_WIND)
ERA_INTERIM_1979_2016[:,5]=np.array(M_ssrd_mon[0:-15])
ERA_INTERIM_1979_2016[:,6]=np.array(M_strd_mon[0:-15])
ERA_INTERIM_1979_2016[:,7]=np.array(M_sp)

np.savetxt('ERA_INTERIM_1979_2016_Malawi.txt',ERA_INTERIM_1979_2016,'%5.2f',delimiter='\t')
