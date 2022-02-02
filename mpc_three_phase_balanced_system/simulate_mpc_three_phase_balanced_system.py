"""
# MPC Simulator for Critical (Priority-weighted) Load Restoration Problem
    # Three-phase Balanced Distribution Grid Test System 
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

## Import the Julia models
    # Data preprocessing module
    # load restoration (lr) model module
    
jl.include("data_handler_threephase.jl")
jl.include("lr_threephase_opf_mpc.jl")

## Paramters

current_dir = os.getcwd()
Δt = 5/60
SOC0 = 90
Zbase = 1
Vbase = 4160
Sbase = (Vbase**2)/Zbase
Cbase = 800

## Get the the distribution grid topology

BUSES, LINES, GENERATORS, WINDTURBINES, PVS, STORAGES = jl.load_case_data(datafile="13buscase")

## Wind and PV power forecast data reader class

exo_data_path = os.chdir(os.path.dirname(current_dir)+'/dataset/exogenous_data')

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

Start_index = 31*24*12 + 2*24*12 + 144    # 2019-08-03 12:00  
num_time_steps = 72                       # 6 hours long with 5 minutes interval
dates_object = ForecastReader()
dates = dates_object.available_dates()
dateslist = list()
for time in dates:
    dateslist.append(time)
datetimes_6h = dateslist[Start_index:Start_index+num_time_steps]

## Get the Wind and PV power forecasts

wind_mult = 150       # kW
pv_mult = 300         # kW
num_time_steps_in_a_day = 288
Pwind_forecast_all = np.empty((0, num_time_steps_in_a_day))
Psolar_forecast_all = np.empty((0, num_time_steps_in_a_day))

for datetime in datetimes_6h:    
    forecast_object = ForecastReader()
    forecast = forecast_object.get_forecast(datetime)
    Pwind = wind_mult*np.array(forecast['wind_gen'])
    Pwind_forecast_all = np.append(Pwind_forecast_all, [Pwind], axis=0)
    Psolar = pv_mult*np.array(forecast['pv_gen'])
    Psolar_forecast_all = np.append(Psolar_forecast_all, [Psolar], axis=0)

## Get the load demands

d_P = np.zeros((len(BUSES), num_time_steps_in_a_day))
d_Q = np.zeros((len(BUSES), num_time_steps_in_a_day))
for b in jl.keys(BUSES):
    for t in range(num_time_steps_in_a_day):
        if b == 2:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 3:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q 
        if b == 5:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 6:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 7:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 9:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 10:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 11:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q
        if b == 12:
            d_P[b,t] = BUSES[b].d_P
            d_Q[b,t] = BUSES[b].d_Q

Pd = d_P       # active_power_demand
Qd = d_Q;      # reactive_power_demand

# Total active power demand at each time step

P_d_tot = np.zeros((num_time_steps_in_a_day))
for jdx in range(num_time_steps_in_a_day):
    P_d_tot[jdx] = np.sum(Pd[:,jdx]) 

Q_d_tot = np.zeros((num_time_steps_in_a_day))
for jdx in range(num_time_steps_in_a_day):
    Q_d_tot[jdx] = np.sum(Qd[:,jdx]) 

# Number of buses

num_buses = len(BUSES)

# Number of load buses

load_buses = jl.findall(Pd[:,0] != 0.0)
num_load_buses = len(load_buses)

## Run the MPC

# Parameters

es_soc = 90
mt_energy = 1000
active_power_demanded = Pd
reactive_power_demanded = Qd
active_power_restored = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
reactive_power_restored = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
num_buses = len(BUSES)
num_lines = len(LINES)
num_lines_singlePhase = 3*num_lines
control_horizon = num_time_steps

# Arrays to store results

Pmtt = ([])
Qmtt = ([])
Pess = ([])
Qess = ([])
SOCes = ([])
Pwt_cutt = ([])
Ppv_cutt = ([])
Pwt_instant = ([])
Ppv_instant = ([])
Pd_instant = np.empty((0, num_buses))
Qd_instant = np.empty((0, num_buses))
P_restoredd = np.empty((0, num_buses))
Q_restoredd = np.empty((0, num_buses))
mu_PP = np.empty((0, num_buses))
mu_QQ = np.empty((0, num_buses))
Pwt_instant = ([])
Ppv_instant = ([])
P_line = np.empty((0, num_lines))
Q_line = np.empty((0, num_lines))
volts = np.empty((0, num_buses))

for i in range(num_time_steps):
    
    wt_power = 1000*Pwind_forecast_all[i, 0:72]/Sbase
    pv_power = 1000*Psolar_forecast_all[i, 0:72]/Sbase
    
    objective_value, P_restored, Q_restored, Pmt, Qmt, Pwtb, Pwt_cut, Ppvs,\
    Ppv_cut, Pes, Qes, SOC_es, voltages, mu_P, mu_Q, frombus, tobus, P_lineflow, Q_lineflow =\
    jl.opf_mpc(BUSES, LINES, GENERATORS, WINDTURBINES, PVS, STORAGES, control_horizon, es_soc, mt_energy, wt_power, pv_power, active_power_demanded, reactive_power_demanded, active_power_restored, reactive_power_restored)
    
    # Apply the first control actions and discard the rest
    
    P_restoredd = np.append(P_restoredd, [P_restored[:,0]], axis=0)
    Q_restoredd = np.append(Q_restoredd, [Q_restored[:,0]], axis=0)
    
    mu_PP = np.append(mu_PP, [mu_P[:,0]], axis=0)
    mu_QQ = np.append(mu_QQ, [mu_Q[:,0]], axis=0)
    
    Pmtt = np.append(Pmtt, Pmt[0]) 
    Qmtt = np.append(Qmtt, Qmt[0])
    
    Pess = np.append(Pess, Pes[0])
    Qess = np.append(Qess, Qes[0])
    
    SOCes = np.append(SOCes, SOC_es[0])
    
    Pwt_cutt = np.append(Pwt_cutt, Pwt_cut[0])
                             
    Ppv_cutt = np.append(Ppv_cutt, Ppv_cut[0])
        
    P_line = np.append(P_line, [P_lineflow[:,0]], axis=0)
    Q_line = np.append(Q_line, [Q_lineflow[:,0]], axis=0)
         
    volts = np.append(volts, [voltages[:,0]], axis=0)
    
    # Save the intant (current) values of the loads and RE generations 
        
    Pd_instant = np.append(Pd_instant, [active_power_demanded[:,0]], axis=0)
    Qd_instant = np.append(Qd_instant, [reactive_power_demanded[:,0]], axis=0)
    Pwt_instant = np.append(Pwt_instant, Pwtb[0])
    Ppv_instant = np.append(Ppv_instant, Ppvs[0])
    
    # Update the storage (ES) battery SOC and microturbine (MT) status for the next step
    
    es_soc = SOC_es[0]
    mt_energy = mt_energy - Pmt[0]*Δt
    
    # Update the restored load for the next step
    
    active_power_restored = np.transpose(P_restored)[0]*1000/Sbase
    reactive_power_restored = np.transpose(Q_restored)[0]*1000/Sbase
    
    # Update the control horizon for the next step
    
    control_horizon -= 1
    
## ------------------ End of the MPC loop --------------------------

## Process the MPC results

Pd_instant = np.transpose(Pd_instant)
Qd_instant = np.transpose(Qd_instant)
P_restoredd = np.transpose(P_restoredd)
Q_restoredd = np.transpose(Q_restoredd)
volts = np.transpose(volts)
P_line = np.transpose(P_line)
Q_line = np.transpose(Q_line)
mu_PP = np.transpose(mu_PP)
mu_QQ = np.transpose(mu_QQ)

# Demanded power at each time step

P_demanded_t = ([])
for pdx in range(Pd_instant.shape[1]):
    P_demanded_t.append(sum(Pd_instant[:,pdx])*Sbase/1000)
    
Q_demanded_t = ([])
for qdx in range(Qd_instant.shape[1]):
    Q_demanded_t.append(sum(Qd_instant[:,qdx])*Sbase/1000)

# Restored load at each time step

P_restored_t = ([])
Q_restored_t = ([])
for jdx in range(num_time_steps):
    P_restored_t = np.append(P_restored_t, sum(P_restoredd[:,jdx]))
    Q_restored_t = np.append(Q_restored_t, sum(Q_restoredd[:,jdx]))
    
# Restored energy at each node

E_restored_n = ([])
for qdx in range(num_buses):
    E_restored_n = np.append(E_restored_n, sum(P_restoredd[qdx,:])*Δt)
    
## Print some of the results

print("Restored Loads:", P_restoredd)
print("Aggregate Restored Load:", P_restored_t)
print("Wind Power Dispatch:", Pwt_instant)
print("Solar Power Dispatch:", Ppv_instant)
print("Microturbine Power Dispatch:", Pmtt)
print("Microturbine Reactive Power Dispatch:", Qmtt)
print("Battery Power Dispatch:", Pess)
print("Battery Reactive Power Dispatch:", Qess)
print("Renewable Power Curtailed:", Pwt_cutt + Ppv_cutt)

## Plot some of the results

# Plot 1: DER dispatch and aggregate load restored

plt.figure(figsize=(12, 6))
plt.plot(P_restored_t, label='Restored Load', linewidth = 3)
plt.plot(Pwt_instant, label='Wind', linewidth = 3)
plt.plot(Ppv_instant, label='PV', linewidth = 3)
plt.plot(Pmtt, label='Microturbine', linewidth = 3)
plt.plot(Pess, label='Battery', linewidth = 3, color = 'blueviolet')
plt.plot((Pwt_cutt + Ppv_cutt), label='Curtailed RE', linewidth = 2, color = 'black')
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Power [kW]', fontsize=15)
plt.title('DER Dispatching & Load Restoration', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 2: DER dispatch and aggregate load restored (reactive power)

plt.figure(figsize=(12, 6))
plt.plot(Q_restored_t, label='Restored Load', linewidth = 3)
plt.plot(Qmtt, label='Microturbine', linewidth = 3)
plt.plot(Qess, label='Battery', linewidth = 3, color = 'blueviolet')
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Reactive Power [kvar]', fontsize=15)
plt.title('DER Dispatching & Load Restoration - Reactive Power', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 3: Change in Energy of the battery and microturbine

Emt_gen = ([])
Emt_gen = np.append(Emt_gen, 1000 - Pmtt[0]*Δt)
for v in range(1,num_time_steps):
    Emt_gen = np.append(Emt_gen, Emt_gen[v-1] - Pmtt[v]*Δt)

Cbat = 800    # kWh

plt.figure(figsize=(12, 6))
plt.plot(Emt_gen, label='Microturbine', linewidth = 3)
plt.plot(Cbat*SOCes/100, label='Battery', linewidth = 3)
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Energy [kWh]', fontsize=15)
plt.title('Remaining Energy From Microturbine and Battery', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(fontsize=14)
plt.show()

# Plot 4: Percentage of the restored load at each time step and node
    # Restored Load = (Restored Load / Demanded Load)*100%
    
plt.figure(figsize=(12, 6))
P_restored_t_perc = (P_restoredd / (Pd*Sbase/1000)[:,0:num_time_steps]) *100
plt.plot(P_restored_t_perc[2], label = 'Load1', linewidth = 3)
plt.plot(P_restored_t_perc[3], label = 'Load2', linewidth = 3)
plt.plot(P_restored_t_perc[5], label = 'Load3', linewidth = 3)
plt.plot(P_restored_t_perc[6], label = 'Load4', linewidth = 3)
plt.plot(P_restored_t_perc[7], label = 'Load5', linewidth = 3)
plt.plot(P_restored_t_perc[9], label = 'Load6', linewidth = 3)
plt.plot(P_restored_t_perc[10], label = 'Load7', linewidth = 3)
plt.plot(P_restored_t_perc[11], label = 'Load8', linewidth = 3)
plt.plot(P_restored_t_perc[12], label = 'Load9', linewidth = 3)    
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Restored Loads [%]', fontsize=15)
plt.tick_params(labelsize=12)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.title('Restored Loads in %')
plt.grid()
plt.legend(loc='upper center', bbox_to_anchor=(0.2, 0.3), fancybox=True, shadow=True, ncol=3)
plt.legend(fontsize=14)
plt.show()

# Plot 5: Nodal voltages

lower_vol_limit = 0.95
upper_vol_limit = 1.05
lb_volts = ([])
ub_volts = ([])
for j in range(num_time_steps):
    lb_volts = np.append(lb_volts, lower_vol_limit)
    ub_volts = np.append(ub_volts, upper_vol_limit)

Buses = ['650','646','645','632','633','634','611','684','671','692','675','652','680']

plt.figure(figsize=(12,6))
for bx in range(len(volts)):
    plt.plot(volts[bx], label = Buses[bx], linewidth = 3)
plt.plot(lb_volts,  label = "LB", linewidth = 4, color="red") 
plt.plot(ub_volts, label = "UB", linewidth = 4, color="red") 
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Voltages [pu]', fontsize=15)
plt.tick_params(labelsize=15)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.ylim((0.94, 1.06))
plt.grid()
plt.legend(loc='upper center', bbox_to_anchor=(0.71, 0.85), fancybox=True, shadow=True, ncol=5)
plt.legend(fontsize=14)
plt.show()