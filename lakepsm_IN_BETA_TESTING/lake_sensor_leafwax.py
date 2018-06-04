# Sylvia Dee <sylvia_dee@brown.edu>
# PRYSM
# PSM d18O Leaf Waxes
# SENSOR MODEL
# Function 'wax_sensor'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
#====================================================================

def wax_sensor(dDp,fC_3,fC_4,eps_c3=-112.8,eps_c4=-124.5):
'''
	SENSOR SUB-MODEL for GDGT proxy data

	INPUTS:
	dDp: delta-D of precipitation at location, if known. 
	NOTE: leaf wax sensor model requires isotope-enabled model output.

	fC3 and fC4: the fractions of C3 and C4 plants, respectively, 


	and εalkane−acid is the average alkane acid offset observed between C28 n-acid and 
	C27 n-alkane in n-alkane biosynthesis [Sachse et al., 2012; Konecky et al., 2016].

	Biological fractionation factors are subsampled from a distribution built based on 
	observations reported in Sachse et al. [2012], where εC27.alk−PC3 = −112.8 ± 34.7
	and εC27.alk−PC 4 = −124.5 ± 28.2

	OUTPUT: pseudoproxy (leaf wax) timeseries

'''
	eps = ((fC_3 * eps_c3) +  (fC_4 * eps_c4))
	delta_d_wax = (dDp + 1000) *(eps /1000 +1.0) - 1000
	return pseudoproxy

