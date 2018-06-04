from pylab import *
import sys, os, cdtime, cdutil, cdms2,  MV2, time, datetime
import numpy as np
import matplotlib.pyplot as plt
from math import exp

#import Scientific.IO.NetCDF 

#======================================================================

# Load Datasets.

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

cdutil.setTimeBoundsMonthly(T_T2M)
T_T2M_ac=cdutil.ANNUALCYCLE.climatology(T_T2M)

cdutil.setTimeBoundsMonthly(M_T2M)
M_T2M_ac=cdutil.ANNUALCYCLE.climatology(M_T2M)

#2m Dew#======================================================================

d2m_M= f('d2m',longitude=(33, 36), latitude = (-14, -9))
d2m_T= f('d2m',longitude=(28, 31), latitude = (-8, -3))

d2m_M.shape
d2m_T.shape

T_D2M = cdutil.averager(d2m_T, axis='xy')
M_D2M = cdutil.averager(d2m_M, axis='xy')

cdutil.setTimeBoundsMonthly(T_D2M)
T_D2M_ac=cdutil.ANNUALCYCLE.climatology(T_D2M)

cdutil.setTimeBoundsMonthly(M_D2M)
M_D2M_ac=cdutil.ANNUALCYCLE.climatology(M_D2M)

# U-wind#======================================================================

u10_M= f('u10',longitude=(33, 36), latitude = (-14, -9))
u10_T= f('u10',longitude=(28, 31), latitude = (-8, -3))

u10_M.shape
u10_T.shape

T_u10 = cdutil.averager(u10_T, axis='xy')
M_u10 = cdutil.averager(u10_M, axis='xy')


cdutil.setTimeBoundsMonthly(T_u10)
T_u10_ac=cdutil.ANNUALCYCLE.climatology(T_u10)

cdutil.setTimeBoundsMonthly(M_u10)
M_u10_ac=cdutil.ANNUALCYCLE.climatology(M_u10)

# V-wind#======================================================================

v10_M= f('v10',longitude=(33, 36), latitude = (-14, -9))
v10_T= f('v10',longitude=(28, 31), latitude = (-8, -3))

v10_M.shape
v10_T.shape

T_v10 = cdutil.averager(v10_T, axis='xy')
M_v10 = cdutil.averager(v10_M, axis='xy')


cdutil.setTimeBoundsMonthly(T_v10)
T_v10_ac=cdutil.ANNUALCYCLE.climatology(T_v10)

cdutil.setTimeBoundsMonthly(M_v10)
M_v10_ac=cdutil.ANNUALCYCLE.climatology(M_v10)

# WindSpeed#======================================================================


T_WIND_ac=np.sqrt(np.power(T_u10_ac,2)+np.power(T_v10_ac,2))
M_WIND_ac=np.sqrt(np.power(M_u10_ac,2)+np.power(M_v10_ac,2))

#Surface Pressure#======================================================================

sp_M= f('sp',longitude=(33, 36), latitude = (-14, -9))
sp_T= f('sp',longitude=(28, 31), latitude = (-8, -3))

sp_M.shape
sp_T.shape

T_sp = cdutil.averager(sp_T, axis='xy')
M_sp = cdutil.averager(sp_M, axis='xy')


cdutil.setTimeBoundsMonthly(T_sp)
T_sp_ac=cdutil.ANNUALCYCLE.climatology(T_sp)

cdutil.setTimeBoundsMonthly(M_sp)
M_sp_ac=cdutil.ANNUALCYCLE.climatology(M_sp)

#Wind Speed 10m#======================================================================

si10_M= f('si10',longitude=(33, 36), latitude = (-14, -9))
si10_T= f('si10',longitude=(28, 31), latitude = (-8, -3))

si10_M.shape
si10_T.shape

T_si10 = cdutil.averager(si10_T, axis='xy')
M_si10 = cdutil.averager(si10_M, axis='xy')

cdutil.setTimeBoundsMonthly(T_si10)
T_si10_ac=cdutil.ANNUALCYCLE.climatology(T_si10)

cdutil.setTimeBoundsMonthly(M_si10)
M_si10_ac=cdutil.ANNUALCYCLE.climatology(M_si10)

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

cdutil.setTimeBoundsMonthly(T_ssrd)
T_ssrd_ac=cdutil.ANNUALCYCLE.climatology(T_ssrd)

cdutil.setTimeBoundsMonthly(M_ssrd)
M_ssrd_ac=cdutil.ANNUALCYCLE.climatology(M_ssrd)

#Radiation#======================================================================

strd_M= f('strd',longitude=(33, 36), latitude = (-14, -9))
strd_T= f('strd',longitude=(28, 31), latitude = (-8, -3))

strd_M.shape
strd_T.shape

T_strd = cdutil.averager(strd_T, axis='xy')
M_strd = cdutil.averager(strd_M, axis='xy')

cdutil.setTimeBoundsMonthly(T_strd)
T_strd_ac=cdutil.ANNUALCYCLE.climatology(T_strd)

cdutil.setTimeBoundsMonthly(M_strd)
M_strd_ac=cdutil.ANNUALCYCLE.climatology(M_strd)


#CONVERSIONS======================================================================
# Convert all variables to /month
# 1/s * 1min/60s * 1hr/60min * 1day/24 hours * 1 month*30 days

T_strd_mon_ac=T_strd_ac*(1./60.*1./60.*1./24.)#*1./30.)
M_strd_mon_ac=M_strd_ac*(1./60.*1./60.*1./24.)#*1./30.)
T_ssrd_mon_ac=T_ssrd_ac*(1./60.*1./60.*1./24.)#*1./30.)
M_ssrd_mon_ac=M_ssrd_ac*(1./60.*1./60.*1./24.)#*1./30.)

# Convert Dewpoint temperature to RH

T_RH_ac = 100.-5.*(T_T2M_ac-T_D2M_ac)
M_RH_ac = 100.-5.*(M_T2M_ac-M_D2M_ac)


#Time Axis#======================================================================

days=np.zeros(12)
days[0]=15.
for i in range(len(days)-1):
	days[i+1]=days[i]+np.floor(365./12.)


year_array=np.ones(12)

#Save#======================================================================

ERA_INTERIM_climatology=np.zeros((12,8))
ERA_INTERIM_climatology[:,0]=year_array
ERA_INTERIM_climatology[:,1]=days
ERA_INTERIM_climatology[:,2]=np.array(T_T2M_ac)
ERA_INTERIM_climatology[:,3]=np.array(T_RH_ac)
ERA_INTERIM_climatology[:,4]=np.array(T_WIND_ac)
ERA_INTERIM_climatology[:,5]=np.array(T_ssrd_mon_ac)
ERA_INTERIM_climatology[:,6]=np.array(T_strd_mon_ac)
ERA_INTERIM_climatology[:,7]=np.array(T_sp_ac)

np.savetxt('ERA_INTERIM_climatology_Tanganyika.txt',ERA_INTERIM_climatology,'%5.2f',delimiter='\t')

#Save#======================================================================

ERA_INTERIM_climatology=np.zeros((12,8))
ERA_INTERIM_climatology[:,0]=year_array
ERA_INTERIM_climatology[:,1]=days
ERA_INTERIM_climatology[:,2]=np.array(M_T2M_ac)
ERA_INTERIM_climatology[:,3]=np.array(M_RH_ac)
ERA_INTERIM_climatology[:,4]=np.array(M_WIND_ac)
ERA_INTERIM_climatology[:,5]=np.array(M_ssrd_mon_ac)
ERA_INTERIM_climatology[:,6]=np.array(M_strd_mon_ac)
ERA_INTERIM_climatology[:,7]=np.array(M_sp_ac)

np.savetxt('ERA_INTERIM_climatology_Malawi.txt',ERA_INTERIM_climatology,'%5.2f',delimiter='\t')


