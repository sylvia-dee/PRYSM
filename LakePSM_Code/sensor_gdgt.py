# Sylvia Dee <sylvia_dee@brown.edu>
# PRYSM
# PSM for Lacustrine Sedimentary Archives
# SENSOR MODEL: GDGT-based measurements, e.g. TEX86, MBT5e
# Function 'gdgt_sensor'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
#====================================================================
def gdgt_sensor(LST,MAAT,beta=(1/50),model='TEX86-loomis'):
	'''
	 SENSOR SUB-MODEL for GDGT proxy data

	 INPUTS: 
	 LST:    LAKE SURFACE TEMPERATURE (C)
	 MAAT:   MEAN ANNUAL AIR TEMPERATURE (C)
	 beta:   User-Specified transfer fucntion of LST to TEX/GDGT
	 model:  Published calibrations for GDGTs. Options:
	 			- Tierney et al., 2008 = 'TEX86-tierney'
	 			- Loomis et al., 2012 = 'TEX86-loomis'
	 			- Russell et al., 2018 = 'MBT' (for MBT-5Me calibrations with air temperature)

	 More recent calibration studies have put forward univariate calibrations for 
	 brGDGTs in lakes with mean annual air temperatures, specifically for MBT 
	 indices which use only 5-methyl isomers MBT'_5ME.
      
     OUTPUT: pseudoproxy timeseries (monthly)
	'''

	if model=='beta':
		pseudoproxy= beta*LST
	elif model=='TEX86-tierney':
		pseudoproxy= (LST+3.4992)/38.874 #tierney2008northern
	elif model=='TEX86-loomis':
		pseudoproxy= (LST+10.92)/54.88 #loomis2012calibration
	elif model=='MBT':
		pseudoproxy=(MAAT + 1.21)/32.42 #russell2018distributions

	return pseudoproxy