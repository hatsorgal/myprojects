# -*- coding: utf-8 -*-
"""
Created on Mon May 15 23:57:21 2023

@author: ghats
"""

# Libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from random import gauss

# Functions
def calculate_theta(P_initial):
    return (13.6e6)/(beta*c*P_initial) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )

def calculate_final_angle(results_mean, radiation_length):
    return np.arctan(results_mean[-1] / radiation_length) + Theta_0

def generate_paths(Theta_0, results):
    for ii in range(0, results.shape[0]):
        for jj in range(1, results.shape[1]):
            scatering_theta = gauss(0 , Theta_0)
            delta_yP= np.tan(scatering_theta) * x
            results[ii][jj] = results[ii][jj-1] +delta_yP
    return results

def plot_scatter_to_radiation_length(radiation_length, results_mean, label, degree=100):
    coeffs = np.polyfit(radiation_length, results_mean, degree)
    y_values_poly = np.polyval(coeffs, radiation_length)
    plt.figure()
    plt.plot(radiation_length, results_mean, '.', label = label, markersize = 1.1)
    plt.plot(radiation_length, y_values_poly, '-', label='Trendline')
    plt.xlabel("$Radiation \  Length [m]$")
    plt.ylabel("Scatter")
    plt.legend()
    plt.grid()
    
def plot_final_angles(results):
    final_angles = np.arctan(results[:,-1] / radiation_length[-1]) + Theta_0
    plt.figure()
    plt.hist(final_angles, bins=50)

    plt.xlabel("Final Scattering Angle")
    plt.ylabel("Frequency")
    plt.grid()
    plt.show()

def plot_particle_trajectories(results):
    plt.figure()
    for i in range(results.shape[0]):
        plt.plot(results[i], radiation_length, linewidth=0.3)
    plt.xlabel("Position along x-axis")
    plt.ylabel("Radiation Length")
    plt.title("Particle Trajectories")
    plt.show()



# Initialization
beta = 1
c = 1 
z = 1 
x = 1.1e-3 
x_0=3.54e-3 
x_X_0 = x/x_0
p_initial = 100e9 
Theta_0 = calculate_theta(p_initial)
results = np.zeros((1000,1000))

# Monte Carlo
radiation_length = np.linspace(0, x_X_0*999, 1000)
results = generate_paths(Theta_0, results)
results_mean = np.mean(results, axis=0)

# Plot
plot_scatter_to_radiation_length(radiation_length, results_mean, 'scatter to radiation length')

# Evaluate P from thehta for different initial energies
energies = [50e9, 100e6, 200e6]

for energy in energies:
    Theta_0_energy = calculate_theta(energy)
    results_energy = np.zeros((1000,1000))
    results_energy = generate_paths(Theta_0_energy, results_energy)
    results_mean_energy =np.mean(results_energy, axis=0)
    label = f'scatter to radiation length for {energy/1e6}M[eV] path'
    plot_scatter_to_radiation_length(radiation_length*x_0, results_mean_energy, label)

# Check the Final angle by brown to simulated
Final_Angle = calculate_final_angle(results_mean, radiation_length[-1])
plot_final_angles(results)
plot_particle_trajectories(results)

print(f"Final angle calculated from simulation: {Final_Angle}")