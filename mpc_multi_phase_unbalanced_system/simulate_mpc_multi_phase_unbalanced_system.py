"""
# MPC Simulator for Critical (Priority-weighted) Load Restoration Problem
    # Multi-phase Unbalanced Distribution Grid Test System 
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from datetime import datetime
import os
import math
import cmath

import julia
from julia.api import Julia

jl = Julia(compiled_modules=False)  
jl = julia.Julia()

## Import the Julia models
    # Data preprocessing module
    # load restoration (lr) model module
    
jl.include("data_handler_multiphase.jl")
jl.include("lr_multiphase_opf_mpc.jl")

## Paramters

current_dir = os.getcwd()
Δt = 5/60
SOC0 = 90
Zbase = 1
Vbase = 2400
Sbase = (Vbase**2)/Zbase
Cbase = 800

## Get the the distribution grid topology

BUSES, LINES_SINGLE, LINES_MULTI, GENERATORS, WINDTURBINES, PVS, STORAGES = jl.load_case_data(datafile="13buscase")

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
        if b == 7-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 8-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 9-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 10-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 11-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 12-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q     
        if b == 19-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 20-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 21-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 23-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 25-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 26-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q 
        if b == 27-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 28-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 29-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 30-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 31-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 34-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q
        if b == 35-1:
            d_P[b,t] = BUSES[b+1].d_P
            d_Q[b,t] = BUSES[b+1].d_Q

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

es_soc = [90, 90, 90]
mt_energy = 1000
active_power_demanded = Pd
reactive_power_demanded = Qd
active_power_restored = np.array([0] * num_time_steps)
reactive_power_restored = np.array([0] * num_time_steps)
num_buses = len(BUSES)
num_lines = len(LINES_SINGLE)
control_horizon = num_time_steps

# Arrays to store MPC results

Pmtt = ([])
Qmtt = ([])
Pess1 = ([])
Pess2 = ([])
Pess3 = ([])
Qess1 = ([])
Qess2 = ([])
Qess3 = ([])
SOCes1 = ([])
SOCes2 = ([])
SOCes3 = ([])
Pwt_instant = ([])
Pwt_cutt = ([])
Qwt_invv1 = ([])
Qwt_invv2 = ([])
Qwt_invv3 = ([])
Ppv_instant = ([])
Ppv_cutt = ([])
Qpv_invv1 = ([])
Qpv_invv2 = ([])
Qpv_invv3 = ([])
Pd_instant = np.empty((0, num_buses))
Qd_instant = np.empty((0, num_buses))
P_restoredd = np.empty((0, num_buses))
Q_restoredd = np.empty((0, num_buses))
mu_PP = np.empty((0, num_buses))
mu_QQ = np.empty((0, num_buses))
Pwt_instant = ([])
Ppv_instant = ([])
from_node = np.empty((num_lines))
to_node = np.empty((num_lines))
P_line = np.empty((0, num_lines))
Q_line = np.empty((0, num_lines))
volts = np.empty((0, num_buses))

for i in range(num_time_steps):
        
    wt_power = 1000*Pwind_forecast_all[i, 0:72]/Sbase
    pv_power = 1000*Psolar_forecast_all[i, 0:72]/Sbase
    
    objective_value, P_restored, Q_restored, Pmt, Qmt, Pwtb, Pwt_cut,\
    Qwt_inv_ph1, Qwt_inv_ph2, Qwt_inv_ph3, Ppvs, Ppv_cut, Qpv_inv_ph1,\
    Qpv_inv_ph2, Qpv_inv_ph3, Pes1, Pes2, Pes3, Qes1, Qes2, Qes3,\
    SOC_es1, SOC_es2, SOC_es3, voltages, mu_P, mu_Q, fromnode, tonode, Plineflow, Qlineflow =\
    jl.opf_mpc(BUSES, LINES_SINGLE, LINES_MULTI, GENERATORS, WINDTURBINES, PVS, STORAGES, control_horizon, es_soc, mt_energy,\
                 wt_power, pv_power, active_power_demanded, reactive_power_demanded,
                 active_power_restored, reactive_power_restored)
  
    # Apply the first control actions and discard the rest
    
    P_restoredd = np.append(P_restoredd, [P_restored[:,0]], axis=0)
    Q_restoredd = np.append(Q_restoredd, [Q_restored[:,0]], axis=0)
    
    mu_PP = np.append(mu_PP, [mu_P[:,0]], axis=0)
    mu_QQ = np.append(mu_QQ, [mu_Q[:,0]], axis=0)
    
    Pmtt = np.append(Pmtt, Pmt[0]) 
    Qmtt = np.append(Qmtt, Pmt[0])    
      
    Pess1 = np.append(Pess1, Pes1[0])
    Pess2 = np.append(Pess2, Pes2[0])
    Pess3 = np.append(Pess3, Pes3[0])    
  
    Qess1 = np.append(Qess1, Qes1[0])
    Qess2 = np.append(Qess2, Qes2[0])
    Qess3 = np.append(Qess3, Qes3[0])
        
    SOCes1 = np.append(SOCes1, SOC_es1[0])
    SOCes2 = np.append(SOCes2, SOC_es2[0])
    SOCes3 = np.append(SOCes3, SOC_es3[0])
    
    Pwt_cutt = np.append(Pwt_cutt, Pwt_cut[0])
    Qwt_invv1 = np.append(Qwt_invv1, Qwt_inv_ph1[0])
    Qwt_invv2 = np.append(Qwt_invv2, Qwt_inv_ph2[0])
    Qwt_invv3 = np.append(Qwt_invv3, Qwt_inv_ph3[0])
                         
    Ppv_cutt = np.append(Ppv_cutt, Ppv_cut[0])
    Qpv_invv1 = np.append(Qpv_invv1, Qpv_inv_ph1[0])
    Qpv_invv2 = np.append(Qpv_invv2, Qpv_inv_ph2[0])
    Qpv_invv3 = np.append(Qpv_invv3, Qpv_inv_ph3[0])
  
    P_line = np.append(P_line, [Plineflow[:,0]], axis=0)
    Q_line = np.append(Q_line, [Qlineflow[:,0]], axis=0)
        
    volts = np.append(volts, [voltages[:,0]], axis=0)
    
    # Save the intant (current) values of the loads and RE generations 
        
    Pd_instant = np.append(Pd_instant, [active_power_demanded[:,0]], axis=0)
    Qd_instant = np.append(Qd_instant, [reactive_power_demanded[:,0]], axis=0)
    Pwt_instant = np.append(Pwt_instant, Pwtb[0])
    Ppv_instant = np.append(Ppv_instant, Ppvs[0])
       
    # Update the energy storage (ES) battery SOC and microturbine (MT) fuel status for the next step
    
    es_soc = [SOC_es1[0], SOC_es2[0], SOC_es3[0]]
    mt_energy = mt_energy - Pmt[0]*Δt
    
    # Update the restored load for the next step
    
    active_power_restored = np.transpose(P_restored)[0]*1000/Sbase
    reactive_power_restored = np.transpose(Q_restored)[0]*1000/Sbase
    
    # Update the control horizon for the next step
    
    control_horizon -=1
        
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
print("Microturbine Reactive Power Dispatch - Phase A:", Qess1)
print("Battery Power Dispatch - Phase A::", Pess1)
print("Battery Reactive Power Dispatch - Phase A:", Qess1)
print("Renewable Power Curtailed:", Pwt_cutt + Ppv_cutt)
  
## Plot some of the results

# Plot 1: DER dispatch and aggregate load restored

plt.figure(figsize=(12,6))
plt.plot(P_restored_t, label='Restored Demand', linewidth = 3)
plt.plot(Pwt_instant, label='WT', linewidth = 3)
plt.plot(Ppv_instant, label='PV', linewidth = 3)
plt.plot(Pmtt, label='MT', linewidth = 3)
plt.plot(Pess1, label='ES1', linewidth = 3)
plt.plot(Pess2, label='ES2', linewidth = 3)
plt.plot(Pess3, label='ES3', linewidth = 3)
plt.plot((Pwt_cutt + Ppv_cutt), label='Curtailed Renewable', linewidth = 2, color = 'black')
plt.xlabel('Time [h]', fontsize=15)
plt.ylabel('Power [kW]', fontsize=15)
plt.title('Active Power Dispatch and Load Restoration', fontsize=15)
plt.tick_params(labelsize=14)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(loc='upper center', bbox_to_anchor=(1.13, 0.99), fancybox=True, shadow=True, ncol=1)
plt.legend(fontsize=14)
plt.show()

# Plot 2: Nodal voltages

nodes = ['650.1','650.2','650.3','632.1','632.2','632.3','670.1','670.2','670.3','671.1','671.2',\
         '671.3','680.1','680.2','680.3','633.1','633.2','633.3','634.1','634.2','634.3','645.3',\
         '645.2','646.3','646.2','692.1','692.2','692.3','675.1','675.2','675.3','684.1','684.3',\
         '611.3','652.1']

lower_vol_limit = 0.95
upper_vol_limit = 1.05
lb_volts = ([])
ub_volts = ([])
for j in range(num_time_steps):
    lb_volts = np.append(lb_volts, lower_vol_limit)
    ub_volts = np.append(ub_volts, upper_vol_limit)
  
plt.figure(figsize=(12,6))
for bx in range(len(volts)):    
    plt.plot(volts[bx], label = nodes[bx], linewidth = 3)  
plt.plot(lb_volts,  label = "LB", linewidth = 4, color="red") 
plt.plot(ub_volts, label = "UB", linewidth = 4, color="red") 
plt.xlabel('Time', fontsize=15)
plt.ylabel('Voltages [pu]', fontsize=15)
plt.tick_params(labelsize=14)
plt.xticks(range(0, 73, 12), [str(int(c/12)) + ":00" for c in range(12*12, 18*12+1, 12)])
plt.rcParams["font.family"] = "Times New Roman"
plt.grid()
plt.legend(loc='upper center', bbox_to_anchor=(0.71, 0.85), fancybox=True, shadow=True, ncol=5)
plt.legend(fontsize=14)
plt.show()