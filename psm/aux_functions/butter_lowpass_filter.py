# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 14:23:45 2014

@author: jeg
"""
from scipy.signal import butter, lfilter, filtfilt
import numpy as np

def butter_lowpass(cutoff, fs, order=3):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def filter(x, cutoff, fs=1.0, order=3):
    b, a = butter_lowpass(cutoff, fs, order=order)
    
    m = np.mean(x)
    y = filtfilt(b, a, x-m) + m
    return y