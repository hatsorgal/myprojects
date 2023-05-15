# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:38:51 2023

@author: ghats
"""
import numpy as np
import math
from random import gauss
from random import uniform



#-- Need to recheck it and run for 50 100 200 MeV


class SimulationMaker():
    def __init(beta , c , z , x , x_0 , p_initial ):
        x_X_0 = x/x_0
        Theta_0 = (13.6e6)/(beta*c*p_initial) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
        results = np.zeros((1000,1000))
        
        
        
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
            