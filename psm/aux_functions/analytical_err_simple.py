# Sylvia Dee <sdee@usc.edu>
# PRYSM 
# OBSERVATION SUB-MODEL
# Function 'ANALYTICAL ERROR (SIMPLE VERSION)'
# Modified [sdee] 1/29/2015

def analytical_err_simple(X,sigma=0.1):
    """
    Function 'analytical_err_simple' produces a d18O-cellulose record which takes into account measurement precision and errors.
       
       Input Arguments:
           1. data vector field, final output of sensor and archive model
           2. sigma (assumed precision of measurement/instrumental error). Default is 0.1 permil.
           NOTE: to include your own measurement precision, simply put it in the sigma field.

       Output: 
        Returns two Numpy Arrays, ['upper', 'lower'], which are saved in the main program, 
        and which now include analytic error upper and lower bounds due to measurement precision.
      
        Example of use: Call function 'analytical_err_simple' to put error envelope around dcell data.

        >>upper, lower = analytical_err_simple(dripwater, 0.1) 
    """
    upper = X+sigma
    lower = X-sigma
    return (upper, lower)