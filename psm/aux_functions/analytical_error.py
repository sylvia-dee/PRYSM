# Sylvia Dee <sdee@usc.edu>
# PRYSM 
# OBSERVATION SUB-MODEL
# Function 'ANALYTICAL ERROR'
# Modified [sdee] 11/11/2014
# Cleaned up [jeg] 11/17/2014

#======================================================================================

def analytical_error(X,sigma=0.1,nsamples=10):
    '''
    Function 'analytical_error' adds Gaussian white noise to simulate
    analytical / measurement error.
       
       Input Arguments:
           1. X, numpy array (1D)
           2. sigma (assumed precision of measurement/instrumental error). Default is 0.1.
           NOTE: to include your own measurement precision or change units, simply modify sigma.
           3. nsamples: number of noise realizations (default = 10)

    Output: 
        Returns a Numpy array, Xn which is saved in the main program, 
        and which now includes analytic error due to measurement precision.
         
         Example of use:  simulate error envelope around speleothem proxy data.
         proxy_record = analytical_error(dripwater, 0.1)
    '''
# Add Gaussian Noise
    import numpy as np
    
    noise = np.random.normal(loc=0.0, scale=sigma, size=(X.shape[0],nsamples))
    Xn    = np.tile(np.expand_dims(X, axis=1),(1,nsamples)) + noise
    
# OR 2-sigma interval (TBD)


    return Xn