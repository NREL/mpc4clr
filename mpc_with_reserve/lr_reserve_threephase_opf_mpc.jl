# Priority-weighted Network-aware Load Restoration Problem for Improved Resilience in Distribution Systems

# Install/add relevant Julia packages

#using Pkg
#Pkg.add("JuMP")
#Pkg.add("Cbc")
#Pkg.add("PyCall")

# Import/load relevant Julia packages

using JuMP
using Cbc
using PyCall

# opf-mpc = Network-aware (OPF-driven) MPC for Distribution Grid Optimal Control/Dispatch

function opf_mpc(buses, lines, generators, windturbines, pvs, storages, control_horizon, es_soc, mt_energy, wt_power, pv_power, active_power_demanded, reactive_power_demanded, active_power_restored, reactive_power_restored)

    # Rearrange the inputs to the opf

    Pdemand = active_power_demanded
    Prestore = active_power_restored
    Qdemand = reactive_power_demanded
    Qrestore = reactive_power_restored
    num_time_steps = control_horizon

    # Parameters

    num_buses = length(buses)
    load_buses = findall(Pdemand[:,1] .!= 0.0)
    num_load_buses = length(load_buses)
    Δt = 5/60
    rho_mt = 1.0  # expected probability of the MT  reserve (ramping) product being used in the distribution system
    rho_es = 1.0  # expected probability of the ES reserve (ramping) product being used in the distribution system
    c = 75        # % system ramping requirement as percentage of renewable forecast

    load_priority = [0.0, 0.0, 1.0, 1.0, 0.0, 0.9, 0.85, 0.8, 0.0, 0.65, 0.45, 0.4, 0.3]
    alpha =  0.2          #$/kWh -- penalty for wind power curtailment
    beta = 0.2            #$/kWh -- penalty for PV power curtailment
    psi = 100             #$/kWh -- penalty (to relaxed hard constraints) for violation to the objective function
    phi = 50

    # Base values

    Zbase = 1
    Vbase = 4160
    Sbase = (Vbase^2)/Zbase
    Cbase = 800

    # Index of system components

    bus_set = collect(keys(buses))
    line_set = collect(keys(lines))
    gen_set = collect(keys(generators))
    wind_set = collect(keys(windturbines))
    pv_set = collect(keys(pvs))
    stor_set = collect(keys(storages))

    # Control steps

    time_set = []
    for t in 1:num_time_steps
        time_set = push!(time_set, t)
    end

    # Bus where each component is connected to

    gen_bus = []
    for g in keys(generators)
        push!(gen_bus, generators[g].bus_idx)
    end

    wind_bus = []
    for w in keys(windturbines)
        push!(wind_bus, windturbines[w].bus_idx)
    end

    pv_bus = []
    for p in keys(pvs)
        push!(pv_bus, pvs[p].bus_idx)
    end

    stor_bus = []
    for s in keys(storages)
        push!(stor_bus, storages[s].bus_idx)
    end

    # Root bus - through which the distribution grid is connected with the upstream grid

    root_bus = []
    for b in keys(buses)
        if buses[b].is_root
            push!(root_bus, b)
        end
    end

    # Lines each node coupled to

    lines_to = Dict()
    for l in keys(lines)
        lines_to[lines[l].to_node] = lines[l]
    end

    num_lines = length(lines_to)

    # Voltage (in per unit) at the root bus

    v_root = 1

    # Create an empty/abstract optimization model

    m = Model()
    set_optimizer(m, Cbc.Optimizer)

    # Define the decision variables

    @variables m begin
        v[bus_set, time_set]           # variable for voltage square
        fp[bus_set, time_set]          # variable for active power flow
        fq[bus_set, time_set]          # variable for reactive power flow

        gp[bus_set, time_set]          # variable for active power generation
        gpr[bus_set, time_set]         # variable for ramping up reserve product (active power)
        gq[bus_set, time_set]          # variable for reactive power generation

        spc[bus_set, time_set]         # variable for energy storage (ES) active power charging
        sqc[bus_set, time_set]         # variable for ES reactive power absorption
        spd[bus_set, time_set]         # variable for ES active power discharging
        sqd[bus_set, time_set]         # variable for ES reactive power injection
        ssoc[bus_set, time_set]        # variable for ES state of charge
        spdr[bus_set, time_set]        # variable for ramping up reserve product (active power)

        wdp[bus_set, time_set]         # variable for active power control (Wind conveter)
        #wdq[bus_set, time_set]         # variable for reactive injection/absorption (Wind converter)
        wdpc[bus_set, time_set]        # variable for Wind power curtailment

        pvp[bus_set, time_set]         # variable for active power control (PV Inverter)
        #pvq[bus_set, time_set]         # variable for reactive injection/absorption (PV Inverter)
        pvpc[bus_set, time_set]        # variable for PV power curtailment

        rlp[bus_set, time_set]         # variable for restored load (kW)
        rlq[bus_set, time_set]         # variable for restored load (kvar)
        mup[bus_set, time_set]         # variable for hard constraint relaxation
        muq[bus_set, time_set]         # variable for hard constraint relaxation

        r[time_set]                    # variable for System-wide ramp-up reserve capability requirements
        tau[time_set]                  # variable for hard constraint relaxation for the reserve
    end

    @variable(m, sc[bus_set, time_set], Bin)     # variable for ES charging indicator and inverter reactive power consumption
    @variable(m, sd[bus_set, time_set], Bin)     # variable for ES discharging indicator and inverter reactive power injection indicator

    # Objective function
        # maximize the total /priority-weighted/ load (active and reactive) served
        # penalize restored load shedding
        # penalize wind and pv power curtailment

    @objective(m, Max, Δt*sum(load_priority[b] .* rlp[b,t] for b in load_buses for t in time_set)
                     + Δt*sum(load_priority[b] .* rlq[b,t] for b in load_buses for t in time_set)
                     - Δt*psi*sum(load_priority[b] .* mup[b,t] for b in load_buses for t in time_set)
                     - Δt*psi*sum(load_priority[b] .* muq[b,t] for b in load_buses for t in time_set)
                     - Δt*phi*sum(tau[t] for t in time_set)
                     - Δt*alpha*sum(wdpc[b,t] for b in wind_bus for t in time_set)
                     - Δt*beta*sum(pvpc[b,t] for b in pv_bus for t in time_set))

    # Constraints

    for b in gen_bus
        for t in time_set
    # All buses with generators
            @constraint(m, gp[b,t] <= buses[b].generator.g_P_max)  #Upper limit constraint for active power generation
            @constraint(m, gpr[b,t] <= buses[b].generator.g_P_max)
            @constraint(m, gp[b,t] >= 0)                           #Lower limit constraint for active power generation
            @constraint(m, gpr[b,t] >= 0)                                   #Lower limit constraint for reserve product
            @constraint(m, gq[b,t] <= buses[b].generator.g_Q_max)  #Upper limit constraint for reactive power generation
            @constraint(m, gq[b,t] >= 0)                           #Lower limit constraint for reactive power generation
        end
    end

    for b in gen_bus
        for t in time_set
            @constraint(m, gp[b,t]+gpr[b,t] <= buses[gen_bus[1]].generator.g_P_max)
        end
    end

    for b in gen_bus
    # All buses with generators - fuel usage /energy production/ limit
        @constraint(m, (sum(gp[b,t] for t in time_set) + sum(rho_mt*gpr[b,t] for t in time_set))*Sbase*Δt/1000 <= mt_energy)
        @constraint(m, sum(gq[b,t] for t in time_set)*Sbase*Δt/1000 <= 0.75*mt_energy)
    end

    for b in wind_bus
        for t in time_set
    # All buses with wind turbines
            @constraint(m, wdp[b,t] == wt_power[t])
            #@constraint(m, wdq[b,t] <= sqrt((buses[b].wind.w_S_max)^2 - (wt_power[t])^2))
            #@constraint(m, wdq[b,t] >= -sqrt((buses[b].wind.w_S_max)^2 - (wt_power[t])^2))
            @constraint(m, wdpc[b,t] <= wt_power[t])
            @constraint(m, wdpc[b,t] >= 0)
        end
    end

    for b in pv_bus
        for t in time_set
    # All buses with PVs
            @constraint(m, pvp[b,t] == pv_power[t])
            #@constraint(m, pvq[b,t] <= sqrt((buses[b].pv.p_S_max)^2 - (pv_power[t])^2))
            #@constraint(m, pvq[b,t] >= -sqrt((buses[b].pv.p_S_max)^2 - (pv_power[t])^2))
            @constraint(m, pvpc[b,t] <= pv_power[t])
            @constraint(m, pvpc[b,t] >= 0)
        end
    end

    for b in stor_bus
        for t in time_set
    # All buses with storages
            @constraint(m, spc[b,t] <= sc[b,t]*buses[b].storage.s_P_max)       #Upper limit constraint for charging power
            @constraint(m, spc[b,t] >= 0)                                       #Lower limit constraint for charging power
            @constraint(m, sqc[b,t] <= sc[b,t]*buses[b].storage.s_Q_max)
            @constraint(m, sqc[b,t] >= 0)

            @constraint(m, spd[b,t] <= sd[b,t]*buses[b].storage.s_P_max)
            @constraint(m, spdr[b,t] <= sd[b,t]*buses[b].storage.s_P_max)
            @constraint(m, spd[b,t]+spdr[b,t] <= sd[b,t]*buses[b].storage.s_P_max)
            @constraint(m, spd[b,t] >= 0)
            @constraint(m, spdr[b,t] >= 0)
            @constraint(m, sqd[b,t] <= sd[b,t]*buses[b].storage.s_Q_max)
            @constraint(m, sqd[b,t] >= 0)

            @constraint(m, sc[b,t] + sd[b,t] <= 1)

            @constraint(m, ssoc[b,t] <= buses[b].storage.s_SOC_max/100)
            @constraint(m, ssoc[b,t] >= buses[b].storage.s_SOC_min/100)

            @constraint(m, ssoc[b,t] - (rho_es*spdr[b,t])/(buses[b].storage.s_eff_dischar/100) >= buses[b].storage.s_SOC_min/100)

            if t == 1
                @constraint(m, ssoc[b,t] == es_soc/100 + ((((buses[b].storage.s_eff_char/100)*spc[b,t])/buses[b].storage.s_cap) -
                (spd[b,t]/((buses[b].storage.s_eff_dischar/100)*buses[b].storage.s_cap))))
            end
            if t > 1
                @constraint(m, ssoc[b,t] == ssoc[b,t-1] + (((buses[b].storage.s_eff_char/100)*spc[b,t]/buses[b].storage.s_cap) -
                (spd[b,t]/((buses[b].storage.s_eff_dischar/100)*buses[b].storage.s_cap))))
            end
        end
    end

    for b in load_buses
        for t in time_set
    # All load buses - loads picked up
            @constraint(m, rlp[b,t] >= 0)
            @constraint(m, rlp[b,t] <= Pdemand[b,t])
            @constraint(m, rlq[b,t] >= 0)
            @constraint(m, rlq[b,t] <= Qdemand[b,t])
            # power factor constraint: Qres/Pres = Qdem/Pdem
            @constraint(m, rlq[b,t]*Pdemand[b,t] == rlp[b,t]*Qdemand[b,t])
        end
    end

    for b in load_buses
        for t in time_set
    # All load buses - relaxation
            @constraint(m, mup[b,t] >= 0)
            @constraint(m, muq[b,t] >= 0)
            if t == 1
                @constraint(m, mup[b,t] >= Prestore[b] - rlp[b,t])
                @constraint(m, muq[b,t] >= Qrestore[b] - rlq[b,t])
            end
            if t>1
                @constraint(m, mup[b,t] >= rlp[b,t-1] - rlp[b,t])
                @constraint(m, muq[b,t] >= rlq[b,t-1] - rlq[b,t])
            end
        end
    end

    for b in setdiff(bus_set, gen_bus)
        for t in time_set
    # All buses without generator
            @constraint(m, gp[b,t] == 0)
            @constraint(m, gpr[b,t] == 0)
            @constraint(m, gq[b,t] == 0)
        end
    end

    for b in setdiff(bus_set, wind_bus)
        for t in time_set
    # All buses without wind
            @constraint(m, wdp[b,t] == 0)
            #@constraint(m, wdq[b,t] == 0)
            @constraint(m, wdpc[b,t] == 0)
        end
    end

    for b in setdiff(bus_set, pv_bus)
        for t in time_set
    # All buses without pv
            @constraint(m, pvp[b,t] == 0)
            #@constraint(m, pvq[b,t] == 0)
            @constraint(m, pvpc[b,t] == 0)
        end
    end

    for b in setdiff(bus_set, stor_bus)
        for t in time_set
    # All buses without storage
            @constraint(m, spc[b,t] == 0)
            @constraint(m, sqc[b,t] == 0)
            @constraint(m, spd[b,t] == 0)
            @constraint(m, spdr[b,t] == 0)
            @constraint(m, sqd[b,t] == 0)
            @constraint(m, ssoc[b,t] == 0)
            @constraint(m, sc[b,t] == 0)
            #@constraint(m, scq[b,t] == 0)
            @constraint(m, sd[b,t] == 0)
            #@constraint(m, sdq[b,t] == 0)
        end
    end

    for b in setdiff(bus_set, load_buses)
        for t in time_set
    # All load buses
            @constraint(m, rlp[b,t] == 0)
            @constraint(m, rlq[b,t] == 0)
            @constraint(m, mup[b,t] == 0)
            @constraint(m, muq[b,t] == 0)
        end
    end

    for b in bus_set
        for t in time_set
    # All buses
            @constraint(m, rlp[b,t] - gp[b,t] - (wdp[b,t]-wdpc[b,t]) - (pvp[b,t]-pvpc[b,t]) + spc[b,t] - spd[b,t] + sum(fp[k,t] for k in buses[b].children) == fp[b,t])
            #@constraint(m, rlq[b,t] - gq[b,t] - wdq[b,t] - pvq[b,t] + sqc[b,t] - sqd[b,t] + sum(fq[k,t] for k in buses[b].children) == fq[b,t])
            @constraint(m, rlq[b,t] - gq[b,t] + sqc[b,t] - sqd[b,t] + sum(fq[k,t] for k in buses[b].children) == fq[b,t])

            @constraint(m, v[b,t] <= buses[b].v_max)
            @constraint(m, v[b,t] >= buses[b].v_min)
        end
    end

    for b in setdiff(bus_set, root_bus)
        for t in time_set
    # All buses without root
            b_ancestor = buses[b].ancestor[1]
            @constraint(m, v[b,t] == v[b_ancestor,t] - 2*(lines_to[b].r * fp[b,t] + lines_to[b].x * fq[b,t]))
            #@constraint(m, lines_to[b].s_max^2 >= fp[b]^2 + fq[b]^2)
        end
    end

    # System-wide ramp-up reserve capability requirement
    for t in time_set
        if t == 1
            @constraint(m, r[t] == 0.0)
        end
        if t > 1
            @constraint(m, r[t] == (c/100)*(sum(wdp[b,t] for b in wind_bus)+sum(pvp[b,t] for b in pv_bus)))
        end
    end

    # Penalize dropping reserves
    for t in time_set
        @constraint(m, tau[t] >= 0)
        #if t == 1
        #    tau[t] == 0
        #end
        #if t>1
        @constraint(m, tau[t] >= r[t]-(sum(gpr[b,t] for b in gen_bus)+sum(spdr[b,t] for b in stor_bus)))
        #@constraint(m, tau[t] == r[t]-(sum(gpr[b,t] for b in gen_bus)+sum(spdr[b,t] for b in stor_bus)))
    end

    #Voltage and flow constraints for root node(s)
    for b in root_bus
        for t in time_set
            @constraint(m, v[b,t] == v_root)
            @constraint(m, fp[b,t] == 0)
            @constraint(m, fq[b,t] == 0)
        end
    end

    # Solve the model

    optimize!(m)

    # Extract and prepare the results/solution

    objective_value = getobjectivevalue(m)

    bus_results_P = DataFrame(bus = Any[], time = Any[], rl_P = Any[], g_P = Any[], g_Pr = Any[],
    wind_P = Any[], wind_Pcut = Any[], pv_P = Any[], pv_Pcut = Any[], s_Pch = Any[], s_Pdch = Any[], s_Pdchr = Any[], s_SOC = Any[], mu_PP = Any[])

    bus_results_Q = DataFrame(bus = Any[], time = Any[], rl_Q = Any[], g_Q = Any[], s_Qch = Any[], s_Qdch = Any[], mu_QQ = Any[])
    #bus_results_Q = DataFrame(bus = Any[], time = Any[], rl_Q = Any[], g_Q = Any[], wind_Q = Any[], pv_Q = Any[], s_Qch = Any[], s_Qdch = Any[])

    bus_results_V = DataFrame(bus = Any[], time = Any[], voltage = Any[])

    line_results = DataFrame(line = Any[], from = Any[], to = Any[], time = Any[], f_P = Any[], f_Q = Any[], current = Any[])

    reserve_required = DataFrame(time = Any[], reserve_req = Any[])

    for b in bus_set
        for t in time_set
            resP = [b, t, getvalue(rlp[b,t])*Sbase/1000, getvalue(gp[b,t])*Sbase/1000, getvalue(gpr[b,t])*Sbase/1000, getvalue(wdp[b,t])*Sbase/1000,
            getvalue(wdpc[b,t])*Sbase/1000, getvalue(pvp[b,t])*Sbase/1000, getvalue(pvpc[b,t])*Sbase/1000, getvalue(spc[b,t])*Sbase/1000,
            getvalue(spd[b,t])*Sbase/1000, getvalue(spdr[b,t])*Sbase/1000, getvalue(ssoc[b,t])*100, getvalue(mup[b,t])*Sbase/1000]
            push!(bus_results_P, resP)
            resQ = [b, t, getvalue(rlq[b,t])*Sbase/1000, getvalue(gq[b,t])*Sbase/1000, getvalue(sqc[b,t])*Sbase/1000,
            getvalue(sqd[b,t])*Sbase/1000, getvalue(muq[b,t])*Sbase/1000]
            push!(bus_results_Q, resQ)
            resV = [b, t, sqrt(getvalue(v[b,t]))]
            push!(bus_results_V, resV)
        end
    end

    for t in time_set
        res_req = [t, getvalue(r[t])*Sbase/1000]
        push!(reserve_required, res_req)
    end

    sort!(bus_results_P)
    sort!(bus_results_Q)
    sort!(bus_results_V)
    sort!(line_results)
    sort!(reserve_required)

    P_restored = zeros(num_buses, num_time_steps)
    mu_P = zeros(num_buses, num_time_steps)
    for b in bus_set
        P_restored[b,:] = bus_results_P[!,"rl_P"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
        mu_P[b,:] = bus_results_P[!,"mu_PP"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
    end

    Q_restored = zeros(num_buses, num_time_steps)
    mu_Q = zeros(num_buses, num_time_steps)
    for b in bus_set
        Q_restored[b,:] = bus_results_Q[!,"rl_Q"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
        mu_Q[b,:] = bus_results_Q[!,"mu_QQ"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
    end

    #Pmt_1p = zeros(length(gen_bus), num_time_steps)
    #Qmt_1p = zeros(length(gen_bus), num_time_steps)
    #m = 1
    #for b in gen_bus
    #    Pmt_1p[m,:] = bus_results_P["g_P"][(gen_bus[m]-1)*num_time_steps+1:(gen_bus[m]-1)*num_time_steps+num_time_steps]
    #    Qmt_1p[m,:] = bus_results_Q["g_Q"][(gen_bus[m]-1)*num_time_steps+1:(gen_bus[m]-1)*num_time_steps+num_time_steps]
    #    m +=1
    #end
    #Pmt = zeros(1, num_time_steps)
    #Qmt = zeros(1, num_time_steps)
    #for mdx in 1:num_time_steps
    #    Pmt[mdx] = sum(Pmt_1p[:,mdx])
    #    Qmt[mdx] = sum(Qmt_1p[:,mdx])
    #end

    Pmt = bus_results_P[!,"g_P"][(gen_bus[1]-1)*num_time_steps+1:(gen_bus[1]-1)*num_time_steps+num_time_steps]
    Qmt = bus_results_Q[!,"g_Q"][(gen_bus[1]-1)*num_time_steps+1:(gen_bus[1]-1)*num_time_steps+num_time_steps]
    Pmt_reserve = bus_results_P[!,"g_Pr"][(gen_bus[1]-1)*num_time_steps+1:(gen_bus[1]-1)*num_time_steps+num_time_steps]

    #Pwtb_1p = zeros(length(wind_bus), num_time_steps)
    #Pwt_cut_1p = zeros(length(wind_bus), num_time_steps)
    #w = 1
    #for b in wind_bus
    #    Pwtb_1p[w,:] = bus_results_P["wind_P"][(wind_bus[w]-1)*num_time_steps+1:(wind_bus[w]-1)*num_time_steps+num_time_steps]
    #    Pwt_cut_1p[w,:] = bus_results_P["wind_Pcut"][(wind_bus[w]-1)*num_time_steps+1:(wind_bus[w]-1)*num_time_steps+num_time_steps]
    #    w +=1
    #end
    #Pwtb = zeros(1, num_time_steps)
    #Pwt_cut = zeros(1, num_time_steps)
    #for wdx in 1:num_time_steps
    #    Pwtb[wdx] = sum(Pwtb_1p[:,wdx])
    #    Pwt_cut[wdx] = sum(Pwt_cut_1p[:,wdx])
    #end

    Pwtb = bus_results_P[!,"wind_P"][(wind_bus[1]-1)*num_time_steps+1:(wind_bus[1]-1)*num_time_steps+num_time_steps]
    Pwt_cut = bus_results_P[!,"wind_Pcut"][(wind_bus[1]-1)*num_time_steps+1:(wind_bus[1]-1)*num_time_steps+num_time_steps]


    #Ppvs_1p = zeros(length(pv_bus), num_time_steps)
    #Ppv_cut_1p = zeros(length(pv_bus), num_time_steps)
    #p = 1
    #for b in pv_bus
    #    Ppvs_1p[p,:] = bus_results_P["pv_P"][(pv_bus[p]-1)*num_time_steps+1:(pv_bus[p]-1)*num_time_steps+num_time_steps]
    #    Ppv_cut_1p[p,:] = bus_results_P["pv_Pcut"][(pv_bus[p]-1)*num_time_steps+1:(pv_bus[p]-1)*num_time_steps+num_time_steps]
    #    p +=1
    #end
    #Ppvs = zeros(1, num_time_steps)
    #Ppv_cut = zeros(1, num_time_steps)
    #for pdx in 1:num_time_steps
    #    Ppvs[pdx] = sum(Ppvs_1p[:,pdx])
    #    Ppv_cut[pdx] = sum(Ppv_cut_1p[:,pdx])
    #end

    Ppvs = bus_results_P[!,"pv_P"][(pv_bus[1]-1)*num_time_steps+1:(pv_bus[1]-1)*num_time_steps+num_time_steps]
    Ppv_cut = bus_results_P[!,"pv_Pcut"][(pv_bus[1]-1)*num_time_steps+1:(pv_bus[1]-1)*num_time_steps+num_time_steps]


    #Qwt_inv = bus_results_Q["wind_Q"][(wind_bus[1]-1)*num_time_steps+1:(wind_bus[1]-1)*num_time_steps+num_time_steps]
    #Qpv_inv = bus_results_Q["pv_Q"][(pv_bus[1]-1)*num_time_steps+1:(pv_bus[1]-1)*num_time_steps+num_time_steps]

    Pes_char = -1 .* bus_results_P[!,"s_Pch"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]
    Pes_dischar = bus_results_P[!,"s_Pdch"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]
    Pes = Pes_char + Pes_dischar
    Pes_reserve = bus_results_P[!,"s_Pdchr"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]

    Qes_inv_cons = -1 .* bus_results_Q[!,"s_Qch"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]
    Qes_inv_inj = bus_results_Q[!,"s_Qdch"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]
    Qes = Qes_inv_cons + Qes_inv_inj

    SOC_es = bus_results_P[!,"s_SOC"][(stor_bus[1]-1)*num_time_steps+1:(stor_bus[1]-1)*num_time_steps+num_time_steps]

    reserve_requiredd = reserve_required[!,"reserve_req"]
    reserve_product = Pmt_reserve + Pes_reserve

    voltages = zeros(num_buses, num_time_steps)
    for bdx in 1:num_buses
        voltages[bdx,:] = bus_results_V[!,"voltage"][1+(bdx-1)*num_time_steps: bdx*num_time_steps]
    end

    frombus = zeros(num_lines)
    tobus = zeros(num_lines)
    P_lineflow = zeros(num_lines, num_time_steps)
    Q_lineflow = zeros(num_lines, num_time_steps)
    #for ldx in 1:num_lines
        #frombus[ldx] = line_results["from"][1+(ldx-1)*num_time_steps]
        #tobus[ldx] = line_results["to"][1+(ldx-1)*num_time_steps]
    #    P_lineflow[ldx,:] = line_results["f_P"][1+(ldx-1)*num_time_steps:ldx*num_time_steps]
    #    Q_lineflow[ldx,:] = line_results["f_Q"][1+(ldx-1)*num_time_steps:ldx*num_time_steps]
    #end

    return objective_value, P_restored, Q_restored, Pmt, Qmt, Pwtb, Pwt_cut, Ppvs, Ppv_cut, Pes, Qes, SOC_es, voltages, reserve_requiredd, reserve_product, mu_P, mu_Q
    #return objective_value, P_restored, Q_restored, Pmt, Qmt, Pwtb, Pwt_cut, Qwt_inv, Ppvs, Ppv_cut, Qpv_inv, Pes, Qes, SOC_es, voltages
end
