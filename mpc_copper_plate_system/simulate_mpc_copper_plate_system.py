"""
# MPC Simulator for Critical (Priority-weighted) Load Restoration Problem
    # Copper Plate Distribution Grid Test System 
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os

import julia
from julia.api import Julia

jl = Julia(compiled_modules=False)  
jl = julia.Julia()

## Forecast (Wind & PV) Data Reader 

class ForecastReader:
    def __init__(self, forecast_pointers_file='forecast_pointers.csv'):
        self._forecast_pointers = pd.read_csv(forecast_pointers_file, parse_dates=['Time'],
                                        index_col=0, infer_datetime_format=True).to_dict()['forecast_file']

    def get_forecast(self, datetime):
        '''
        Get the forecast as a pandas dataframe for the datetime
        '''
        k = pd.Timestamp(datetime)

        if k not in self._forecast_pointers:
            raise Exception(f"Datetime {datetime} does not have a forecast available")

        return pd.read_csv(self._forecast_pointers[k], parse_dates=['Time'], index_col=0, infer_datetime_format=True)

    def available_dates(self):
        return self._forecast_pointers.keys()
    
## Get the datetimes over the control horizon

current_dir = os.getcwd()
os.chdir(os.path.dirname(current_dir) + "/dataset/exogenous_data")
Start_index = 31*24*12 + 6*24*12 + 72      # 2019-08-07 06:00
num_time_steps = 72                        # 6 hours long with 5 minutes interval
dates_object = ForecastReader()
dates = dates_object.available_dates()
dateslist = list()
for time in dates:
    dateslist.append(time)
datetimes_6h = dateslist[Start_index:Start_index+num_time_steps]

## Get the PV and wind power forecasts

wind_max = 150      # kW
pv_max = 300        # kW
num_time_steps_in_a_day = 288
Pwind_forecast_all = np.empty((0, num_time_steps_in_a_day))
Psolar_forecast_all = np.empty((0, num_time_steps_in_a_day))

for datetime in datetimes_6h:    
    forecast_object = ForecastReader()
    forecast = forecast_object.get_forecast(datetime)
    Pwind = wind_max*np.array(forecast['wind_gen'])
    Pwind_forecast_all = np.append(Pwind_forecast_all, [Pwind], axis=0)
    Psolar = pv_max*np.array(forecast['pv_gen'])
    Psolar_forecast_all = np.append(Psolar_forecast_all, [Psolar], axis=0)

## Get the load demands
    
Load_dataset = pd.read_csv('Constant_Load_Time_Series.csv', header=0, infer_datetime_format=True, 
                           parse_dates=['Node'], index_col=['Node'])
Load_dataset = np.array(Load_dataset)
Pl_demanded = Load_dataset

## Import the Julia/JuMP load restoration optimization model

os.chdir(current_dir)

jl.include("load_restoration_model.jl") 

## Run the MPC

NUM_OF_LOAD_BUS = 10
num_nodes = NUM_OF_LOAD_BUS
Δt = 5/60
battery_SOC = 0.90
mt_energy = 1000
control_horizon = num_time_steps
Pl_restored = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

P_restoredd = np.empty((0, num_nodes))
Pmt_genn = ([])
Pbatt = ([])
SOCbatt = ([])
Pwt_cutt = ([])
Ppv_cutt = ([])
Pwt_instant = ([])
Ppv_instant = ([])
muuu = ([])
Pl_instant = np.empty((0, num_nodes))

for i in range(num_time_steps):               
    
    # Get the wind and PV power at this step

    Pwt = Pwind_forecast_all[i, 0:72]
    Ppv = Psolar_forecast_all[i, 0:72]

    # Run the optimization model
    
    P_restored, Pmt_gen, Pbat, SOCbat, Pwt_cut, Ppv_cut, Objective_value = \
        jl.execute_mpc(control_horizon, battery_SOC, mt_energy, Pwt, Ppv, Pl_demanded, Pl_restored)
    
    # Apply the first control actions and discard the rest
    
    P_restoredd = np.append(P_restoredd, [P_restored[:,0]], axis=0)
    Pmt_genn = np.append(Pmt_genn, Pmt_gen[0])
    Pbatt = np.append(Pbatt, Pbat[0])
    SOCbatt = np.append(SOCbatt, SOCbat[0])
    Pwt_cutt = np.append(Pwt_cutt, Pwt_cut[0])
    Ppv_cutt = np.append(Ppv_cutt, Ppv_cut[0])    
    
    # Save the intant (current) values of the loads and RE generations 
        
    Pl_instant = np.append(Pl_instant, [Pl_demanded[:,0]], axis=0)
    Pwt_instant = np.append(Pwt_instant, Pwt[0])
    Ppv_instant = np.append(Ppv_instant, Ppv[0])
    
    # Update the energy storage (ES) battery SOC and microturbine (MT) fuel status for the next step
    
    battery_SOC = SOCbat[0]
    mt_energy = mt_energy - Pmt_gen[0]*Δt
    
    # Update the restored load for the next step
    
    Pl_restored = np.transpose(P_restored)[0]
    
    # Update the control horizon for the next step

    control_horizon -= 1

## ------------------ End of the MPC loop --------------------------


## Process the MPC results

# Demanded and restored loads

Pl_instant = np.transpose(Pl_instant)
P_restoredd = np.transpose(P_restoredd)
    
# Demanded load (kW) at each time step 

Plbase = np.array([33, 34, 8.5, 85, 60, 60, 58, 115, 64, 85])
Pldemanded_t = list()
for i in range(len(Plbase)):
    Pldemanded_t.append([Plbase[i]] * num_time_steps)
Pldemanded_t = np.array(Pldemanded_t)

Actual_Load_Real_t = list()
for ix in range(num_time_steps):
    Load_t = sum(Pl_instant[:,ix])
    Actual_Load_Real_t.append(Load_t)
Actual_Load_Real_t = np.array(Actual_Load_Real_t)

# Demanded energy (kWh) at each time step 

Energy_n = list()
for jx in range(len(Pl_instant)):
    ener_n = sum(Pl_instant[jx,:])*Δt
    Energy_n.append(ener_n)
Energy_n = np.array(Energy_n)

# Aggregate restored load at each time step

P_restored_t = list()
for jdx in range(num_time_steps):
    P_restored_t.append(sum(P_restoredd[:,jdx]))
P_restored_t = np.array(P_restored_t)

# Aggregate restored energy at each node

E_restored_n = list()
for qdx in range(num_nodes):
    E_restored_n.append((sum(P_restoredd[qdx,:]))*Δt)
E_restored_n = np.array(E_restored_n)

## Print some of the results

print("Restored Loads:", P_restoredd)
print("Aggregate Restored Load:", P_restored_t)
print("Wind Power Dispatch:", Pwt_instant)
print("Solar Power Dispatch:", Ppv_instant)
print("Microturbine Power Dispatch:", Pmt_genn)
print("Battery Power Dispatch:", Pbatt)
print("Renewable Power Curtailed:", Pwt_cutt + Ppv_cutt)

## Plot some of the results

# Plot 1: DER dispatch and aggregate load restored

plt.figure(figsize=(12, 6))
plt.plot(P_restored_t, label='Restored Load', linewidth = 3)
plt.plot(Pwt_instant, label='Wind', linewidth = 3)
plt.plot( Ppv_instant, label='PV', linewidth = 3)
plt.plot(Pmt_genn, label='Microturbine', linewidth = 3)
plt.plot(Pbatt, label='Battery', linewidth = 3, color = 'blueviolet')
plt.plot((Pwt_cutt + Ppv_cutt), label='Curtailed RE', linewidth = 2, color = 'black')
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Power [kW]', fontsize=15)
plt.title('DER Dispatching & Load Restoration', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(6*12, 12*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 2: Change in Energy of the battery and microturbine

Emt_gen = list()
Emt_gen.append(1000 - Pmt_genn[0]*Δt)
for v in range(1,num_time_steps):
    Emt_gen.append(Emt_gen[v-1] - Pmt_genn[v]*Δt)
Emt_gen = np.array(Emt_gen)

Cbat = 800    # kWh

plt.figure(figsize=(12, 6))
plt.plot(Emt_gen, label='Microturbine', linewidth = 3)
plt.plot(Cbat*SOCbatt, label='Battery', linewidth = 3)
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Energy [kWh]', fontsize=15)
plt.title('Remaining Energy From Microturbine and Battery', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(6*12, 12*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 3: Percentage of the restored load at each time step and node
    # Restored Load = (Restored Load / Demanded Load)*100%
    
plt.figure(figsize=(12, 6))
P_restored_t_perc = (P_restoredd / Pl_demanded[:,0:72]) *100
Nodes = ['Node1','Node2','Node3','Node4','Node5','Node6','Node7','Node8','Node9','Node10']
for dx in range(len(P_restored_t_perc)):
    plt.plot(P_restored_t_perc[dx], label = Nodes[dx], linewidth = 3)
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Restored Loads [%]', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(6*12, 12*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.title('Restored Loads in %')
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 4: Actual/demanded energy vs restored energy per node

plt.figure(figsize=(12, 6))
plt.plot(Energy_n, label='Demanded', linewidth = 3)
plt.plot( E_restored_n, label='Served', linewidth = 3)
plt.plot((Energy_n - E_restored_n), label='Shedded', linewidth = 3)
plt.xlabel('Nodes - Priority of Loads Decreases From Node 1 to 10', fontsize=15)
plt.ylabel('Energy [kWh]', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 10, 1), [str(int(c)) for c in range(1, 11, 1)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()