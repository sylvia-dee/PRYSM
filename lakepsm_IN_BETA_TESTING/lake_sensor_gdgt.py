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
	 			- Powers et al., 2011 = 'TEX86-powers'
	 			- Loomis et al., 2012 = 'TEX86-loomis'
	 			- Russell et al., 2018 = 'MBT-R' (for MBT-5Me calibrations with air temperature)
				- De Jonge et al., 2014 = 'MBT-J'
	 More recent calibration studies have put forward univariate calibrations for 
	 brGDGTs in lakes with mean annual air temperatures, specifically for MBT 
	 indices which use only 5-methyl isomers MBT'_5ME.
      
     OUTPUT: pseudoproxy timeseries (monthly)
	'''

	if model=='beta':
		pseudoproxy= beta*LST
	elif model=='TEX86-tierney':
		pseudoproxy= (LST+3.4992)/38.874 #tierney2008northern
	elif model=='TEX86-powers':
		pseudoproxy=(LST+14.0)/55.0   #ALST = 55.0*TEX86 − 14.0, with an r2 = 0.86, corresponding to an error of the calibration of ± 3.7 °C
	elif model=-'TEX86-loomis':
		pseudoproxy= (LST+10.92)/54.88 #loomis2012calibration
	elif model=='MBT-R':
		pseudoproxy=(MAAT + 1.21)/32.42 #russell2018distributions
	elif model=='MBT-J':
		pseudoproxy=(MAAT + 8.57)/31.45 # de jonge 2014
    return pseudoproxy