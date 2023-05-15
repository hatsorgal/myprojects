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
energies = [50e6, 100e6, 200e6]

for energy in energies:
    Theta_0_energy = calculate_theta(energy)
    results_energy = np.zeros((1000,1000))
    results_energy = generate_paths(Theta_0_energy, results_energy)
    results_mean_energy =np.mean(results_energy, axis=0)
    label = f'scatter to radiation length for {energy/1e6}M[eV] path'
    plot_scatter_to_radiation_length(radiation_length*x_0, results_mean_energy, label)

# Check the Final angle by brown to simulated
Final_Angle = calculate_final_angle(results_mean, radiation_length[-1])
print(f"Final angle calculated from simulation: {Final_Angle}")







#%% evaluate P from thehta

#----- Generate 50 Gev P_inital paths to get-------

P_initial_50 = 50e9
Theta_0_50 = (13.6e6)/(beta*c*P_initial_50) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results_50 = np.zeros((1000,1000))
scatering_theta = np.zeros(1000)
for ii in range(1,len(results_50)):
    for jj in range(1,len(results_50)):
        scatering_theta[ii] = scatering_theta[ii] + gauss(0 , Theta_0_50)
        delta_yP= np.tan(scatering_theta[ii])*x
        results_50[ii][jj]= results_50[ii][jj-1] +delta_yP



results_mean_50 =np.zeros(1000)  
for ii in range(0,len(results)):
    count=np.zeros(1000)  
    for jj in range(1,len(results)):
        count[jj-1] =results_50[ii][jj-1]
    results_mean_50[ii] = np.mean(count)+results_mean_50[ii-1]   
    
    
Final_Angle_50=np.zeros(1000)
P_final_estimated_50 =np.zeros(1000)

for jj in range(1,len(results_50)):
    Final_Angle_50[jj] = math.atan(results_50[len(results_50)-1][jj] / 
                            radiation_length[len(results_50)-1] ) 
    P_final_estimated_50[jj] = (13.6e6)/(beta*c*np.abs(Final_Angle_50[jj])) *z*np.sqrt(1/300) * ( 1 + 0.038*math.log(1/300) ) / 1e9
     
    

plt.figure()
plt.plot(radiation_length*x_0 ,results_mean_50,'+', label = 'scatter to radiation' + 
         'length for 50M[eV] path', markersize=0.3 )
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()




#----- Generate  100  Mev P_inital paths to get-------

P_initial_100 = 100e6
Theta_0_100 = (13.6e6)/(beta*c*P_initial_100) *z*np.sqrt(x_X_0) * ( 1 + 0.038*math.log(x_X_0) )
results_100 = np.zeros((1000,1000))
for ii in range(1,len(results_100)):
    for jj in range(1,len(results_100)):
        scatering_theta = gauss(0 , Theta_0_100)
        #delta_yP = scatering_theta * x  # Small angles -> tanx ~ x
        results_100[ii][jj]=results_100[ii][jj-1] +scatering_theta #+delta_yP

Final_Angle_100=np.zeros(1000)
P_final_estimated_100 =np.zeros(1000)

for jj in range(1,len(results_100)):
    Final_Angle_100[jj] = math.atan(results_100[jj][len(results_100)-1] / 
                            radiation_length[len(results_100)-1] )*180/np.pi
    P_final_estimated_100[jj] = (13.6e6)/((beta*c*np.abs(Final_Angle_100[jj])) *z*np.sqrt(1/300) * ( 1 + 0.038*math.log(1/300) ) 
     *1e6)
     
plt.figure()       
plt.hist(P_final_estimated_100, bins= 10000)  
plt.xlim(0, 3000)
plt.xlabel("Final P estimated")
plt.ylabel("Count")    


plt.figure()
plt.plot(radiation_length ,results_100,'.', label = 'scatter to radiation' + 
         'length for 100M[eV] path', markersize=1.1)
plt.xlabel("$Radiation \  Length [m]$")
plt.ylabel("Scatter")
plt.legend()
plt.grid()

Final_Angle_100 = math.atan(results_100[len(results_100)-1] / 
                        radiation_length[len(results_100)-1] ) 
P_final_estimated_100 = np.abs((13.6e6)/(beta*c*Final_Angle_100) *z*np.sqrt(1/300) * ( 1 + 0.038*math.log(1/300) )/1e6)
 
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
 






    
