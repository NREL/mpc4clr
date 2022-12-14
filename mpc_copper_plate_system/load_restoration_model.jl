# Priority-weighted Load Restoration Problem for Improved Resilience in Distribution Systems

# Install/add relevant Julia packages

using Pkg
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
#Pkg.add("GLPK")
#Pkg.add("Cbc")
#Pkg.add("SCIP")
#Pkg.add("Alpine")
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("OpenDSSDirect")
#Pkg.add("Plots")
#Pkg.add("PyPlot")
#Pkg.add("GR")
#Pkg.add("UnicodePlots")
#Pkg.add("PyCall")

# Import/load relevant Julia packages

using JuMP
using Ipopt
using GLPK
using Cbc
using CSV
using DataFrames
using OpenDSSDirect
using Plots
using PyCall
#using SCIP
#using Alpine
#using OSQP
#using TimeSeries
#using ParameterJuMP
#using Statistics

# execute_mpc() takes the battery available SOC, microturbine (MT) available fuel/energy,
# and forecasts of the load demand and renewables (wind and PV)
# and returns an optimized actions of battery (charge/discharge) power, MT power,
# and wind and pv curtailments.

function execute_mpc(control_horizon, battery_SOC, mt_energy, Pwt, Ppv, Pl_demanded, Pl_restored)

    # Define the constant paramters

    NUM_OF_LOAD_BUS = 10
    num_nodes = NUM_OF_LOAD_BUS
    num_time_steps = control_horizon
    FIVE_MIN_TO_HOUR = 5/60
    Δt = FIVE_MIN_TO_HOUR
    load_priority_weight = [1.0, 1.0, 0.9, 0.85, 0.8, 0.65, 0.45, 0.4, 0.3, 0.3]
    PmtMIN = 0            #kW
    PmtMAX = 300
    PchMAX = 200
    PdischMAX = 200
    SOCMIN = 0.20          #%
    SOCMAX = 1.0
    eff_ch = 0.95
    eff_disch = 0.90
    Ces = 800              #kWh
    alpha =  0.2          #$/kWh -- penalty for wind power curtailment
    beta = 0.2            #$/kWh -- penalty for PV power curtailment
    psi = 100      #$/kWh -- penalty (to relaxed hard constraints) for violation to the objective function

    # Create an empty optimization model

    model = Model()
    set_optimizer(model, Cbc.Optimizer)

    # Define the decision variables

    @variables model begin
        Pl[1:num_nodes, 1:num_time_steps]
        Pmt[1:num_time_steps]
        Pch[1:num_time_steps]
        Pdisch[1:num_time_steps]
        SOC[1:num_time_steps]
        Pwt_curt[1:num_time_steps]
        Ppv_curt[1:num_time_steps]
        mu[1:num_nodes, 1:num_time_steps]
    end

    @variable(model, bch[1:num_time_steps], Bin)
    @variable(model, bdisch[1:num_time_steps], Bin);

    # Add constraints

        # Power balance

    @constraint(model, [j=1:num_time_steps],
                (Pwt[j]-Pwt_curt[j]) + (Ppv[j]-Ppv_curt[j]) + Pmt[j] - Pch[j] + Pdisch[j]
        == sum(Pl[1:num_nodes, j]))

        # Load Picked up

    @constraint(model, [i=1:num_nodes, j=1:num_time_steps],
                Pl[i, j] >= 0)
    @constraint(model, [i=1:num_nodes, j=1:num_time_steps],
                Pl[i, j] <= Pl_demanded[i, j])

    @constraint(model, [i=1:num_nodes, j=1],
                mu[i, j] >= Pl_restored[i] - Pl[i, j])
    # mu[i, j] >= Pl_restored[i] - Pl[i, j]

        # Relaxation

    @constraint(model, [i=1:num_nodes, j=2:num_time_steps],
                mu[i, j] >= Pl[i, j-1] - Pl[i, j])
     # max((Pl_demanded[i, j-1] - Pl_demanded[i, j]), 0)
    @constraint(model, [i=1:num_nodes, j=1:num_time_steps],
                mu[i, j] >= 0)

        # Microturbine

    @constraint(model, [j=1:num_time_steps],
                Pmt[j] >= PmtMIN)
    @constraint(model, [j=1:num_time_steps],
                Pmt[j] <= PmtMAX)
    @constraint(model,(sum(Pmt[1:num_time_steps])) * Δt <= mt_energy)

        # Energy storage

    @constraint(model, [j=1:num_time_steps],
                Pch[j] >= 0)
    @constraint(model, [j=1:num_time_steps],
                Pch[j] <= bch[j] * PchMAX)
    @constraint(model, [j=1:num_time_steps],
                Pdisch[j] >= 0)
    @constraint(model, [j=1:num_time_steps],
                Pdisch[j] <= bdisch[j] * PdischMAX)
    @constraint(model, [j=1:num_time_steps],
                bch[j] + bdisch[j] == 1)
    @constraint(model, [j=1:num_time_steps],
                SOC[j] >= SOCMIN)
    @constraint(model, [j=1:num_time_steps],
                SOC[j] <= SOCMAX)
    @constraint(model, [j=1],
                SOC[j] == battery_SOC + ((eff_ch*Pch[j]/Ces) - (Pdisch[j]/(eff_disch*Ces))) * Δt)
    @constraint(model, [j=2:num_time_steps],
                SOC[j] == SOC[j-1] + ((eff_ch*Pch[j]/Ces) - (Pdisch[j]/(eff_disch*Ces))) * Δt)

     # Curtailment -- renewable energy generation

    @constraint(model, [j=1:num_time_steps],
                Pwt_curt[j] >= 0)
    @constraint(model, [j=1:num_time_steps],
                Pwt_curt[j] <= Pwt[j])
    @constraint(model, [j=1:num_time_steps],
                Ppv_curt[j] >= 0)
    @constraint(model, [j=1:num_time_steps],
                Ppv_curt[j] <= Ppv[j])

    # Objective function: maximize the total priority-weighted loads picked up

     @objective(model, Max,
             Δt*sum(load_priority_weight[1:num_nodes] .* Pl[1:num_nodes, 1:num_time_steps])
           - Δt*psi*sum(load_priority_weight[1:num_nodes] .* mu[1:num_nodes, 1:num_time_steps])
           - Δt*alpha*sum(Pwt_curt[1:num_time_steps])
           - Δt*beta*sum(Ppv_curt[1:num_time_steps]))

# - Δt*psi*sum(load_priority_weight[1:num_nodes] .* mu[1:num_nodes, 1:num_time_steps])

    # Solve

    optimize!(model)

    # Extract the model solution

    P_restored = JuMP.value.(Pl)
    Pmt_gen = JuMP.value.(Pmt)
    Pbat_charge = -1 .* JuMP.value.(Pch)
    Pbat_discharge = JuMP.value.(Pdisch)
    Pbat = Pbat_charge .+ Pbat_discharge
    SOCbat = JuMP.value.(SOC)
    Is_Charging = JuMP.value.(bch)
    Is_Discharging = JuMP.value.(bdisch)
    Pwt_cut = JuMP.value.(Pwt_curt)
    Ppv_cut = JuMP.value.(Ppv_curt)
    muu = JuMP.value.(mu)
    Objective_value = objective_value(model)

    return P_restored, Pmt_gen, Pbat, SOCbat, Pwt_cut, Ppv_cut, Objective_value
end
