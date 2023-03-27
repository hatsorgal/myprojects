# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%% Functions
import math
import numpy as np # math functions
import scipy # scientific functions
import matplotlib.pyplot as plt # for plotting figures and setting their properties
import pandas as pd # handling data structures (loaded from files)
from scipy.stats import linregress # contains linregress (for linear regression)
from scipy.optimize import curve_fit as cfit # non-linear curve fitting
from sklearn.metrics import r2_score # import function that calculates R^2 score
from scipy.signal import resample
from random import gauss
from random import uniform





#%% Initial configuration
beta = 1
c = 1 # eV
z = -1 ### charge of the particle
x_X_0 = uniform(1e-3,1e2) ### radiation length gcm^-2
p_initial = 29.7885e6 # in eV
Theta_0 = (13.6e6)/(beta*c*p_initial) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results = np.ones(100)


#%% Monte Carlo 
temp = Theta_0

for ii in range(0,len(results)):
    results[ii] = temp;
    temp = gauss(-3*temp, 3*temp)
    
    
    
