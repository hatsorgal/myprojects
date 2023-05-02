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
x = 1.1e-3 #[m]
x_0=3.54e-3 #[m]
x_X_0 = x/x_0
p_initial = 100e6 # in eV
Theta_0 = (13.6e6)/(beta*c*p_initial) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results = np.zeros((1000,1000))


#%% Monte Carlo 
temp = Theta_0
radiation_length = np.zeros(1000);
radiation_length[0]=x_X_0
for k in range(1,len(radiation_length)):
    radiation_length[k] = k* x_X_0 ;
    
for ii in range(0,len(results)):
    for jj in range(1,len(results)):
        scatering_theta = gauss(0 , Theta_0)
        delta_yP= np.tan(scatering_theta) * x
        results[ii][jj] = results[ii][jj-1] +delta_yP

results_mean =np.zeros(1000)  
for ii in range(0,len(results)):
    count=np.zeros(1000)  
    for jj in range(1,len(results)):
        count[jj-1] =results[ii][jj-1]
    results_mean[ii] = np.mean(count)+results_mean[ii-1]       
#%% Some graphs, need to ask Enrique what

plt.figure()
plt.plot(radiation_length ,results_mean,'.', label = 'scatter to radiation length')
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()

plt.figure()       
plt.hist(results_mean, bins= 100)  
plt.xlim(-0.004, 0.004)
plt.xlabel("Scatter")
plt.ylabel("Count")


#%% Check the Final angle by brown to simulated

Final_Angle = math.asin(results_mean[len(results_mean)-1] / 
                        radiation_length[len(results_mean)-1] )   + Theta_0



#-------------- Brown motion? which one ? -----------------





#%% evaluate P from thehte

#----- Generate 50 Mev P_inital paths to get-------

P_initial_50 = 50e6
Theta_0_50 = (13.6e6)/(beta*c*P_initial_50) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results_50 = np.zeros((1000,1))
for ii in range(1,len(results_50)):
    scatering_theta = gauss(0 , Theta_0_50)
    delta_yP= np.tan(scatering_theta) * x
    results_50[ii]=results_50[ii-1] +delta_yP

plt.figure()
plt.plot(radiation_length ,results_50,'.', label = 'scatter to radiation' + 
         'length for 50M[eV] path', markersize=1.1)
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()

Final_Angle_50 = math.atan(results_50[len(results_50)-1] / 
                        radiation_length[len(results_50)-1] ) +Theta_0_50
P_final_estimated_50 = (13.6e6)/(beta*c*Final_Angle_50) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
 
#----- Generate  100  Mev P_inital paths to get-------

P_initial_100 = 100e6
Theta_0_100 = (13.6e6)/(beta*c*P_initial_100) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results_100 = np.zeros((1000,1))
for ii in range(1,len(results_100)):
    scatering_theta = gauss(0 , Theta_0_100)
    delta_yP= np.tan(scatering_theta) * x
    results_100[ii]=results_100[ii-1] +delta_yP

plt.figure()
plt.plot(radiation_length ,results_100,'.', label = 'scatter to radiation' + 
         'length for 100M[eV] path', markersize=1.1)
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()

Final_Angle_100 = math.atan(results_100[len(results_100)-1] / 
                        radiation_length[len(results_100)-1] ) +Theta_0_100
P_final_estimated_100 = (13.6e6)/(beta*c*Final_Angle_100) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
 
#----- Generate 200 Mev P_inital paths to get-------

P_initial_200 = 200e6
Theta_0_200 = (13.6e6)/(beta*c*P_initial_200) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results_200 = np.zeros((1000,1))
for ii in range(1,len(results_200)):
    scatering_theta = gauss(0 , Theta_0_200)
    delta_yP= np.tan(scatering_theta) * x
    results_200[ii]=results_200[ii-1] +delta_yP

plt.figure()
plt.plot(radiation_length ,results_200,'.', label = 'scatter to radiation' + 
         'length for 200M[eV] path', markersize=1.1)
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()

Final_Angle_200 = math.atan(results_200[len(results_200)-1] / 
                        radiation_length[len(results_200)-1] ) +Theta_0_200
P_final_estimated_200 = (13.6e6)/(beta*c*Final_Angle_200) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
 






    
