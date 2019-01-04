# Sylvia Dee <sylvia_dee@brown.edu>
# PRYSM
# PSM d18O Leaf Waxes
# SENSOR MODEL
# Function 'wax_sensor'
# Modified 03/8/2016 <sylvia_dee@brown.edu>
#====================================================================

def wax_sensor(dDp,fC_3,fC_4,eps_c3=-112.8,eps_c4=-124.5):

	eps = ((fC_3 * eps_c3) +  (fC_4 * eps_c4))
	delta_d_wax = (dDp + 1000) *(eps /1000 +1.0) - 1000
	return delta_d_wax


def wax_uncertainty(dDp,fC_3,fC_4,eps_c3=-112.8,eps_c4=-124.5,eps_c3_err=34.7, eps_c4_err=124.5):
	C3_MCA=np.linspace(eps_c3-eps_c3_err,eps_c3+eps_c3_err,1000)
	C4_MCA=np.linspace(eps_c4-eps_c4_err,eps_c4+eps_c4_err,1000)

	from numpy import random

	# INITIALIZE MCA ARRAY

	delta_d_wax_mc=np.zeros((1000,np.shape(dDp)))

	# COMPUTE NEW ETOT VALUE 1000 TIMES, SAMPLING UNCERTAINTIES IN EPSILON:

	for k in range(1000):

		E29_c3 = np.random.choice(C3_MCA,1)
		E29_c4 = np.random.choice(C4_MCA,1)
		eps = ((fC_3 * E29_c3) +  (fC_4 * eE29_c4))
		delta_d_wax_mc[i,:] = (dDp + 1000) *(eps /1000 +1.0) - 1000

	Q1=np.percentile(delta_d_wax_mc,2.5,axis=0)
	Q2=np.percentile(delta_d_wax_mc,97.5,axis=0)
	return delta_d_wax_mc,Q1,Q2