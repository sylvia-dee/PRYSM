import pytest

import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))

import numpy as np


# from psm.cellulose import sensor as cs

def test_cellulose():
    pass

     
# t = np.linspace(0,100)
# T = 100
# P = 5
# RH = 80
# d18Os = 1
# d18Op = 0.6
# d18Ov = 0.1


# def test_cellulose_sensor_default_initialize():
    
#     sensor = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov)

#     # here we just want to check that it initializes with the defaults.
#     pass


# def test_cellulose_sensor_run_with_defaults():
    
#     sensor = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov)

#     sensor.run_sensor()

#     assert not np.isnan(sensor.dcell)


# def test_cellulose_sensor_run_with_model_options():
    
#     sensor1 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='roden')
#     sensor2 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='evans')

#     sensor1.run_sensor()
#     sensor2.run_sensor()

#     assert not np.isnan(sensor1.dcell)
#     assert not np.isnan(sensor2.dcell)


# def test_cellulose_sensor_run_model_option_capital_letters():
    

#     sensor1 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='roden')
#     sensor2 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='evans')

#     sensor3 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='Roden')
#     sensor4 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='EvAnS')

#     sensor1.run_sensor()
#     sensor2.run_sensor()
#     sensor3.run_sensor()
#     sensor4.run_sensor()

#     assert sensor1.dcell == sensor3.dcell
#     assert sensor2.dcell == sensor4.dcell

# def test_cellulose_sensor_change_run_roden():
    
#     sensor1 = cs.CelluloseSensor(t, T, P, RH, d18Os, d18Op, d18Ov, model_flag='roden')

#     newval = 100000
#     sensor1.set_RH(newval)

#     assert sensor1._model.RH == newval
#     assert sensor1._model.RH == sensor1._RH

#     sensor1.run_sensor()

#     assert not np.isnan(sensor1.dcell)




# # this is of course not an exhaustive list of tests, you should test all sorts
# # of things, including realistic scenarios, and making sure answers are
# # realistic.
