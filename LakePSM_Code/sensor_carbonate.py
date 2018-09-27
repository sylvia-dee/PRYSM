import numpy

def carb_sensor(LST, d18Ow, isoflag = -1, model = 'ONeil'):
	'''
	 SENSOR SUB-MODEL for Carbonate Proxy Data:Calculation of carbonate isotopic 
	 ratio based on background d18O_w and temp. 

	 INPUTS: 
	 LST:    LAKE SURFACE TEMPERATURE (C)
	 d18Ow:  d18O of lake water

	 model:  Published calibrations for d18O of lacustrine carbonate-temperature relationship. 
	 		Options:
	 			- O'Neil et al., 1969 = 'Oneil' 
	 			- Kim and O'Neil, 1997 = 'Kim-Oneil'
	 			- Erea and Luz, 1983 = 'Erezluz'
	 			- Bemis et al.,1983 = 'Bemis'
	 			- Jean Lynch (need cittation)='Lynch'

     OUTPUT: pseudoproxy timeseries of d18O-carb (monthly)
	'''

	temp = LST
	if isoflag == -1: 
		d18O_w = -2.
	else:
		d18O_w = d18Ow

	# Includes correction to VPDB scale for d18Ow
	if (model == 'ONeil'):
		dw = d18O_w - 0.2
		dOc = dw + 21.9 - numpy.sqrt(21.9*21.9-10.*(16.9-temp))
	elif (model == 'Kim-ONeil'):
		dw = d18O_w - 0.27
		dOc = dw + 25.8 - numpy.sqrt(25.8*25.8 - 11.1*(16.1-temp))
	elif (model == 'ErezLuz'):
		dw = d18O_w - 0.22
		dOc = dw + 75.3 - numpy.sqrt(75.3*75.3 -33.3*(17.0-temp))
	elif (model == 'Bemis'):
		dw = d18O_w - 0.27
		dOc = dw + (13.2-temp)/4.89
	elif (model == 'Lynch'):
		dw = d18O_w - 0.27
		dOc = dw + (16.21 - temp)/4.99

	return dOc