# Multi-phase Distribution Grid

# Priority-weighted Network-aware Load Restoration Problem for Improved Resilience in Distribution Systems

# Install/add relevant Julia packages

#using Pkg
#Pkg.add("JuMP")
#Pkg.add("GLPK")
#Pkg.add("PyCall")

# Import/load relevant Julia packages

using JuMP
using GLPK
using PyCall
using LinearAlgebra
using OpenDSSDirect

# opf-mpc = Network-aware (OPF-driven) MPC for Distribution Grid Optimal Control/Dispatch

function opf_mpc(nodes, lines_singlephase, lines_multiphase, generators, windturbines, pvs, storages, control_horizon, es_soc, mt_energy, wt_power,
                 pv_power, active_power_demanded, reactive_power_demanded, active_power_restored, reactive_power_restored)

    # Rearrange the inputs to the opf

    Pdemand = active_power_demanded
    Prestore = active_power_restored
    Qdemand = reactive_power_demanded
    Qrestore = reactive_power_restored
    num_time_steps = control_horizon

    # Parameters
    #num_bus = 14
    #max_num_phase_per_bus = 3
    num_nodes = length(nodes)
    #load_nodes = findall(Pdemand[:,1] .!= 0.0)
    #num_load_nodes = length(load_nodes)
    Δt = 5/60

    #load_priority = [0.0, 0.0, 1.0, 1.0, 0.9, 0.0, 0.85, 0.8, 0.0, 0.65, 0.45, 0.4, 0.3]
    #load_priority = [0.0, 1.0, 1.0, 1.0, 0.95, 0.95, 0.95, 0.90, 0.0, 0.90, 0.85,
    #                 0.85, 0.85, 0.8, 0.8, 0.8, 0.75, 0.70, 0.6, 0.5, 0.4, 0.0,
    #                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    #load_priority = [1.0, 1.0, 1.0, 0.95, 0.95, 0.95, 0.90, 0.0, 0.90, 0.85,
    #                 0.85, 0.85, 0.8, 0.8, 0.8, 0.75, 0.70, 0.6, 0.5, 0.4, 0.0,
    #                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    load_priority = [0, 0, 0, 0, 0, 0, 1.0, 1.0, 1.0, 0.95, 0.95, 0.95, 0, 0,
                     0, 0, 0, 0, 0.90, 0.85, 0.85, 0, 0.80, 0, 0.75,
                     0.70, 0.70, 0.70, 0.65, 0.60, 0.55, 0, 0, 0.50, 0.45]

    alpha =  0.2          #$/kWh -- penalty for wind power curtailment
    beta = 0.2            #$/kWh -- penalty for PV power curtailment
    psi = 100             #$/kWh -- penalty (to relaxed hard constraints) for violation to the objective function

    # Base values

    Zbase = 1
    Vbase = 2400
    Sbase = (Vbase^2)/Zbase
    Cbase = 800

    ## Admittance Matrix Function

    function form_admittance(Zbase,Vbase, Sbase)

        Nbus = 14
        Ybase = 1/Zbase

        # Impedance matrix

        # Line configurations

        # Configuration 601

        Zs1 = [0.3465+1.0179im 0.1560+0.5017im 0.1580+0.4236im;
               0.1560+0.5017im 0.3375+1.0478im 0.1535+0.3849im;
               0.1580+0.4236im 0.1535+0.3849im 0.3414+1.0348im]/Zbase

        Ys1 = 10^-6*[6.2998im -1.9958im -1.2595im;
                     -1.9958im 5.9597im -0.7417im;
                     -1.2595im -0.7417im 5.6386im]/Ybase

        # Configuration 602

        Zs2 = [0.7526+1.1814im 0.1580+0.4236im 0.1560+0.5017im;
               0.1580+0.4236im 0.7475+1.1983im 0.1535+0.3849im;
               0.1560+0.5017im 0.1535+0.3849im 0.7436+1.2112im]/Zbase

        Ys2 = 10^-6*[5.6990im -1.0817im -1.6905im;
                     -1.0817im 5.1795im -0.6588im;
                     -1.6905im -0.6588im 5.4246im]/Ybase


        # Configuration 603

        Zs3 = [0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
               0.0000+0.0000im 1.3294+1.3471im 0.2066+0.4591im;
               0.0000+0.0000im 0.2066+0.4591im 1.3238+1.3569im]/Zbase

        Ys3 = 10^-6*[0.0000im 0.0000im 0.0000im;
                     0.0000im 4.7097im -0.8999im;
                     0.0000im -0.8999im 4.6658im]/Ybase

        # Configuration 604

        Zs4 = [1.3238+1.3569im 0.0000+0.0000im 0.2066+0.4591im;
               0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
               0.2066+0.4591im 0.0000+0.0000im 1.3294+1.3471im]/Zbase

        Ys4 = 10^-6*[4.6658im 0.0000im -0.8999im;
                     0.0000im 0.0000im 0.0000im;
                     -0.8999im 0.0000im 4.7097im]/Ybase

        # Configuration 605

        Zs5 = [0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
               0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
               0.0000+0.0000im 0.0000+0.0000im 1.3292+1.3475im]/Zbase

        Ys5 = 10^-6*[0.0000im 0.0000im 0.0000im;
                     0.0000im 0.0000im 0.0000im;
                     0.0000im 0.0000im 4.5193im]/Ybase

        # Configuration 606

        Zs6 = [0.7982+0.4463im 0.3192+0.0328im 0.2849-0.0143im;
               0.3192+0.0328im 0.7891+0.4041im 0.3192+0.0328im;
               0.2849-0.0143im 0.3192+0.0328im 0.7982+0.4463im]/Zbase

        Ys6 = 10^-6*[96.8897im 0.0000im 0.0000im;
                     0.0000im 96.8897im 0.0000im;
                     0.0000im 0.0000im 96.8897im]/Ybase

        # Configuration 607

        Zs7 = [1.3425+0.5124im 0.0000+0.0000im 0.0000+0.0000im;
               0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
               0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im]/Zbase

        Ys7 = 10^-6*[88.9912im 0.0000im 0.0000im;
                     0.0000im 0.0000im 0.0000im;
                     0.0000im 0.0000im 0.0000im]/Ybase

        # Line matrices

        convfm = (1/5280)   # 1 mile = 5280 feet

        # Line '650'-'632' linecode = 601
        Z12 = Zs1*(2000)*convfm
        Z12i = pinv(Z12)
        Y12 = 0.5*Ys1*(2000)*convfm

        # Line '632'-'670' lincode = 601
        Z23 = Zs1*(667)*convfm
        Z23i = pinv(Z23)
        Y23 = 0.5*Ys1*(667)*convfm

        # Line '632'-'633' linecode = 602
        Z26 = Zs2*(500)*convfm
        Z26i = pinv(Z26)
        Y26 = 0.5*Ys2*(500)*convfm

        # Line '632'-'645' linecode = 603
        Z28 = Zs3*(500)*convfm
        Z28i = pinv(Z28)
        Y28 = 0.5*Ys3*(500)*convfm

        # Line '670'-'671' linecode = 601
        Z34 = Zs1*(1333)*convfm
        Z34i = pinv(Z34)
        Y34 = 0.5*Ys1*(1333)*convfm

        # Line '671'-'680' linecode = 601
        Z45 = Zs1*(1000)*convfm
        Z45i = pinv(Z45)
        Y45 = 0.5*Ys1*(1000)*convfm

        # Line '671'-'692' linecode = 606
        Z410 = Zs6*(10)*convfm
        Z410i = pinv(Z410)
        Y410 = 0.5*Ys6*(10)*convfm

        # Line '671'-'684' linecode = 604
        Z412 = Zs4*(300)*convfm
        Z412i = pinv(Z412)
        Y412 = 0.5*Ys4*(300)*convfm

        # Line '633'-'634' linecode = 602
        Z67 = Zs2*(10)*convfm
        Z67i = pinv(Z67)
        Y67 = 0.5*Ys2*(10)*convfm

        # Line '645'-'646' linecode = 603
        Z89 = Zs3*(300)*convfm
        Z89i = pinv(Z89)
        Y89 = 0.5*Ys3*(300)*convfm

        # Line '692'-'675' linecode = 606
        Z1011 = Zs6*(500)*convfm
        Z1011i = pinv(Z1011)
        Y1011 = 0.5*Ys6*(500)*convfm

        # Line '684'-'611' linecode = 605
        Z1213 = Zs5*(300)*convfm
        Z1213i = pinv(Z1213)
        Y1213 = 0.5*Ys5*(300)*convfm

        # Line '684'-'652' linecode = 607
        Z1214 = Zs7*(800)*convfm
        Z1214i = pinv(Z1214)
        Y1214 = 0.5*Ys7*(800)*convfm

        # Network admittance matrix

        oo = zeros(3,3)

        Y_net = [Z12i+Y12 -Z12i zeros(3,3*(Nbus-2));
                 -Z12i Z12i+Z23i+Z26i+Z28i+Y12+Y23+Y26+Y28 -Z23i oo oo -Z26i oo -Z28i zeros(3,3*(Nbus-8));
                 oo -Z23i Z23i+Z34i+Y23+Y34 -Z34i zeros(3,3*(Nbus-4));
                 oo oo -Z34i Z34i+Z45i+Z410i+Z412i+Y34+Y45+Y410+Y412 -Z45i oo oo oo oo -Z410i oo -Z412i zeros(3,3*(Nbus-12));
                 oo oo oo -Z45i Z45i+Y45 zeros(3,3*(Nbus-5));
                 oo -Z26i oo oo oo Z26i+Z67i+Y26+Y67 -Z67i zeros(3,3*(Nbus-7));
                 oo oo oo oo oo -Z67i Z67i+Y67 zeros(3,3*(Nbus-7));
                 oo -Z28i oo oo oo oo oo Z28i+Z89i+Y28+Y89 -Z89i zeros(3,3*(Nbus-9));
                 oo oo oo oo oo oo oo -Z89i Z89i+Y89 zeros(3,3*(Nbus-9));
                 oo oo oo -Z410i oo oo oo oo oo Z410i+Z1011i+Y410+Y1011 -Z1011i zeros(3,3*(Nbus-11));
                 oo oo oo oo oo oo oo oo oo -Z1011i Z1011i+Y1011 zeros(3,3*(Nbus-11));
                 oo oo oo -Z412i oo oo oo oo oo oo oo Z412i+Z1213i+Z1214i+Y412+Y1213+Y1214 -Z1213i -Z1214i zeros(3,3*(Nbus-14));
                 oo oo oo oo oo oo oo oo oo oo oo -Z1213i Z1213i+Y1213 zeros(3,3*(Nbus-13));
                 oo oo oo oo oo oo oo oo oo oo oo -Z1214i oo Z1214i+Y1214 zeros(3,3*(Nbus-14))]

        return Y_net
    end

    # Call the form-admittance function
    Ynet = form_admittance(Zbase,Vbase, Sbase)

    Y = Ynet[setdiff(1:end, [22,25,35,37,38,41,42]), setdiff(1:end, [22,25,35,37,38,41,42])]
    YY = Y
    YY[22,:] = Y[23,:]
    YY[23,:] = Y[22,:]
    YY[24,:] = Y[25,:]
    YY[25,:] = Y[24,:]

    # Indeces

    node_set = collect(keys(nodes))

    load_nodes = []
    for ndx in node_set
        if nodes[ndx].d_P != 0.0
            push!(load_nodes, nodes[ndx].index)
        end
    end

    num_nodes = length(nodes)
    num_load_nodes = length(load_nodes)

    bus_set = []
    for i in node_set
        push!(bus_set, nodes[i].bus)
    end
    bus_set_uni = sort(unique(bus_set))

    phase_set = []
    for i in node_set
        push!(phase_set, nodes[i].phase)
    end
    phase_set_uni = sort(unique(phase_set))

    line_set = collect(keys(lines_singlephase))
    gen_set = collect(keys(generators))
    wind_set = collect(keys(windturbines))
    pv_set = collect(keys(pvs))
    stor_set = collect(keys(storages))

    # Control steps

    time_set = []
    for t in 1:num_time_steps
        time_set = push!(time_set, t)
    end

    # Node where each component is connected to

    gen_bus = []
    gen_phase = []
    gen_node = []
    for g in keys(generators)
        push!(gen_bus, generators[g].bus)
        push!(gen_phase, generators[g].phase)
        push!(gen_node, generators[g].node_idx)
    end

    wind_bus = []
    wind_phase = []
    wind_node = []
    for w in keys(windturbines)
        push!(wind_bus, windturbines[w].bus)
        push!(wind_phase, windturbines[w].phase)
        push!(wind_node, windturbines[w].node_idx)
    end

    pv_bus = []
    pv_phase = []
    pv_node = []
    for p in keys(pvs)
        push!(pv_bus, pvs[p].bus)
        push!(pv_phase, pvs[p].phase)
        push!(pv_node, pvs[p].node_idx)
    end

    stor_bus = []
    stor_phase = []
    stor_node = []
    for s in keys(storages)
        push!(stor_bus, storages[s].bus)
        push!(stor_phase, storages[s].phase)
        push!(stor_node, storages[s].node_idx)
    end

    # Root bus - through which the distribution grid is connected with the upstream grid

    root_node = []
    for b in keys(nodes)
        if nodes[b].is_root
            push!(root_node, b)
        end
    end

    # lines each node/bus connected to

    lines_singlephase_to = Dict()
    for l in keys(lines_singlephase)
        lines_singlephase_to[lines_singlephase[l].to_node] = lines_singlephase[l]
    end
    num_lines_singlephase = length(lines_singlephase_to)

    lines_multiphase_to = Dict()
    for l in keys(lines_multiphase)
        lines_multiphase_to[lines_multiphase[l].to_bus] = lines_multiphase[l]
    end
    num_lines_multiphase = length(lines_multiphase_to)


    # Voltage (in per unit) at the root node(s)

    v_root = 1

    # Zbus - System Impedance Matrix

    masterfile = "/Users/aeseye/Desktop/AGM-Resilience-Project/IEEE13bus_for_AGM_Grid_Resiliency/ieee13bus/ieee_13node_DSS/master_islanded3_without_load.dss"
    #masterfile = "/Users/aeseye/Desktop/AGM-Resilience-Project/IEEE13bus_for_AGM_Grid_Resiliency/ieee13bus/ieee_13node_DSS/master_islanded3.dss"
    dss("""
        clear
        compile $masterfile
        solve
    """)

    AllNodeNamesY = Circuit.YNodeOrder()  # Nodes order in the order of SystemY

    Ybus = Circuit.SystemY()
    #Ybus = sparse(Ybus)
    Ybus = Ybus[4:end, 4:end]

    Zbus_self = complex(zeros(length(Ybus[:,1]), length(Ybus[1,:])))
    for row in 1:length(Ybus[:,1])
        for col in 1:length(Ybus[1,:])
            Zbus_self[row,col] = pinv(Ybus[row,col])
            #Zbus_self[row,col] = pinv(YY[row,col])
            #Zbus_self[row,col] = inv(Ybus[row,col])
        end
    end

    #Zbus = inv(Ybus)/Zbase
    #Zbus = pinv(Ybus)/Zbase
    Zbus = Zbus_self/Zbase

    #Zbus = Zm/Zbase
    #Zbus = pinv(YY)
    #Zbus = Zbus_self

    Rb = real(Zbus)
    Xb = imag(Zbus)

    # Create an empty/abstract optimization model

    m = Model()
    set_optimizer(m, GLPK.Optimizer)

    # Define the decision variables

    @variables m begin
        v[node_set, time_set] >= 0      # variable for voltage square

        fp[node_set, time_set]          # variable for active power flow
        fq[node_set, time_set]          # variable for reactive power flow

        gp[node_set, time_set]          # variable for active power generation
        gq[node_set, time_set]          # variable for reactive power generation

        spc[node_set, time_set]         # variable for energy storage (ES) active power charging
        sqc[node_set, time_set]         # variable for ES reactive power absorption
        spd[node_set, time_set]         # variable for ES active power discharging
        sqd[node_set, time_set]         # variable for ES reactive power injection
        ssoc[node_set, time_set]        # variable for ES state of charge

        wdp[node_set, time_set]         # variable for active power control (Wind conveter)
        wdq[node_set, time_set]         # variable for reactive injection/absorption (Wind converter)
        wdpc[node_set, time_set]        # variable for Wind power curtailment

        pvp[node_set, time_set]         # variable for active power control (PV Inverter)
        pvq[node_set, time_set]         # variable for reactive injection/absorption (PV Inverter)
        pvpc[node_set, time_set]        # variable for PV power curtailment

        rlp[node_set, time_set]         # variable for restored load (kW)
        rlq[node_set, time_set]         # variable for restored load (kvar)
        mup[node_set, time_set]         # variable for hard constraint relaxation
        muq[node_set, time_set]         # variable for hard constraint relaxation
    end

    @variable(m, sc[node_set, time_set], Bin)     # variable for ES charging indicator and inverter reactive power consumption
    @variable(m, sd[node_set, time_set], Bin)     # variable for ES discharging indicator and inverter reactive power injection indicator

    # Objective function
        # maximize the total /priority-weighted/ load (active and reactive) served
        # penalize restored load shedding
        # penalize wind and pv power curtailment

    @objective(m, Max, Δt*sum(load_priority[b] .* rlp[b,t] for b in load_nodes for t in time_set)
                     + Δt*sum(load_priority[b] .* rlq[b,t] for b in load_nodes for t in time_set)
                     - Δt*psi*sum(load_priority[b] .* mup[b,t] for b in load_nodes for t in time_set)
                     - Δt*psi*sum(load_priority[b] .* muq[b,t] for b in load_nodes for t in time_set)
                     - Δt*alpha*sum(wdpc[b,t] for b in wind_node for t in time_set)
                     - Δt*beta*sum(pvpc[b,t] for b in pv_node for t in time_set))

    # Constraints

    for b in gen_node
        for t in time_set
    # All nodes with generators
            @constraint(m, gp[b,t] <= nodes[b].generator.g_P_max)  #Upper limit constraint for active power generation
            @constraint(m, gp[b,t] >= 0)                           #Lower limit constraint for active power generation
            @constraint(m, gq[b,t] <= nodes[b].generator.g_Q_max)  #Upper limit constraint for reactive power generation
            @constraint(m, gq[b,t] >= 0)                           #Lower limit constraint for reactive power generation
            #@constraint(m, gq[b,t] >= -nodes[b].generator.g_Q_max)   #Lower limit constraint for reactive power generation
        end
    end

    for b in gen_node
    # All nodes with generators - fuel usage /energy production/ limit
        @constraint(m, sum(gp[b,t] for t in time_set)*Sbase*Δt/1000 <= mt_energy/3)
        @constraint(m, sum(gq[b,t] for t in time_set)*Sbase*Δt/1000 <= 0.75*mt_energy/3)
    end

    for b in wind_node
        for t in time_set
    # All nodes with wind turbines
            @constraint(m, wdp[b,t] == wt_power[t])
            #@constraint(m, wdq[b,t] == 0)
            @constraint(m, wdq[b,t] <= sqrt((nodes[b].wind.w_S_max)^2 - (wt_power[t])^2))
            @constraint(m, wdq[b,t] >= -sqrt((nodes[b].wind.w_S_max)^2 - (wt_power[t])^2))
            @constraint(m, wdpc[b,t] <= wt_power[t])
            @constraint(m, wdpc[b,t] >= 0)
        end
    end

    for b in pv_node
        for t in time_set
    # All nodes with PVs
            @constraint(m, pvp[b,t] == pv_power[t])
            #@constraint(m, pvq[b,t] == 0)
            @constraint(m, pvq[b,t] <= sqrt((nodes[b].pv.p_S_max)^2 - (pv_power[t])^2))
            @constraint(m, pvq[b,t] >= -sqrt((nodes[b].pv.p_S_max)^2 - (pv_power[t])^2))
            @constraint(m, pvpc[b,t] <= pv_power[t])
            @constraint(m, pvpc[b,t] >= 0)
        end
    end

    for b in stor_node
        for t in time_set
    # All nodes with storages
            @constraint(m, spc[b,t] <= sc[b,t]*nodes[b].storage.s_P_max)       #Upper limit constraint for charging power
            @constraint(m, spc[b,t] >= 0)                                      #Lower limit constraint for charging power
            @constraint(m, sqc[b,t] <= sc[b,t]*nodes[b].storage.s_Q_max)
            @constraint(m, sqc[b,t] >= 0)
            #@constraint(m, sqc[b,t] == 0)

            @constraint(m, spd[b,t] <= sd[b,t]*nodes[b].storage.s_P_max)
            @constraint(m, spd[b,t] >= 0)
            @constraint(m, sqd[b,t] <= sd[b,t]*nodes[b].storage.s_Q_max)
            @constraint(m, sqd[b,t] >= 0)
            #@constraint(m, sqd[b,t] == 0)

            #@constraint(m, sc[b,t] + sd[b,t] == 1)
            @constraint(m, sc[b,t] + sd[b,t] <= 1)

            @constraint(m, ssoc[b,t] <= nodes[b].storage.s_SOC_max/100)
            @constraint(m, ssoc[b,t] >= nodes[b].storage.s_SOC_min/100)

            if b == stor_node[1]
                if t == 1
                    @constraint(m, ssoc[b,t] == es_soc[1]/100 + ((((nodes[b].storage.s_eff_char/100)*spc[b,t])/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
                if t > 1
                    @constraint(m, ssoc[b,t] == ssoc[b,t-1] + (((nodes[b].storage.s_eff_char/100)*spc[b,t]/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
            elseif b == stor_node[2]
                if t == 1
                    @constraint(m, ssoc[b,t] == es_soc[2]/100 + ((((nodes[b].storage.s_eff_char/100)*spc[b,t])/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
                if t > 1
                    @constraint(m, ssoc[b,t] == ssoc[b,t-1] + (((nodes[b].storage.s_eff_char/100)*spc[b,t]/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
            elseif b == stor_node[3]
                if t == 1
                    @constraint(m, ssoc[b,t] == es_soc[3]/100 + ((((nodes[b].storage.s_eff_char/100)*spc[b,t])/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
                if t > 1
                    @constraint(m, ssoc[b,t] == ssoc[b,t-1] + (((nodes[b].storage.s_eff_char/100)*spc[b,t]/nodes[b].storage.s_cap)
                    - (spd[b,t]/((nodes[b].storage.s_eff_dischar/100)*nodes[b].storage.s_cap))))
                end
            end
        end
    end

    for b in load_nodes
        for t in time_set
    # All load nodes - loads picked up
            @constraint(m, rlp[b,t] >= 0)
            @constraint(m, rlp[b,t] <= Pdemand[b,t])
            @constraint(m, rlq[b,t] >= 0)
            @constraint(m, rlq[b,t] <= Qdemand[b,t])
            # power factor constraint: Qres/Pres = Qdem/Pdem
            @constraint(m, rlq[b,t]*Pdemand[b,t] == rlp[b,t]*Qdemand[b,t])
        end
    end

    for b in load_nodes
        for t in time_set
    # All load nodes - relaxation
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

    for b in setdiff(node_set, gen_node)
        for t in time_set
    # All nodes without generator
            @constraint(m, gp[b,t] == 0)
            @constraint(m, gq[b,t] == 0)
        end
    end

    for b in setdiff(node_set, wind_node)
        for t in time_set
    # All nodes without wind
            @constraint(m, wdp[b,t] == 0)
            @constraint(m, wdq[b,t] == 0)
            @constraint(m, wdpc[b,t] == 0)
        end
    end

    for b in setdiff(node_set, pv_node)
        for t in time_set
    # All nodes without pv
            @constraint(m, pvp[b,t] == 0)
            @constraint(m, pvq[b,t] == 0)
            @constraint(m, pvpc[b,t] == 0)
        end
    end

    for b in setdiff(node_set, stor_node)
        for t in time_set
    # All nodes without storage
            @constraint(m, spc[b,t] == 0)
            @constraint(m, sqc[b,t] == 0)
            @constraint(m, spd[b,t] == 0)
            @constraint(m, sqd[b,t] == 0)
            @constraint(m, ssoc[b,t] == 0)
            @constraint(m, sc[b,t] == 0)
            @constraint(m, sd[b,t] == 0)
        end
    end

    for b in setdiff(node_set, load_nodes)
        for t in time_set
    # All load nodes
            @constraint(m, rlp[b,t] == 0)
            @constraint(m, rlq[b,t] == 0)
            @constraint(m, mup[b,t] == 0)
            @constraint(m, muq[b,t] == 0)
        end
    end

    for b in node_set
        for t in time_set
    # All nodes
            @constraint(m, rlp[b,t] - gp[b,t] - (wdp[b,t]-wdpc[b,t]) - (pvp[b,t]-pvpc[b,t]) + spc[b,t] - spd[b,t] + sum(fp[k,t] for k in nodes[b].children) == fp[b,t])
            @constraint(m, rlq[b,t] - gq[b,t] - wdq[b,t] - pvq[b,t] + sqc[b,t] - sqd[b,t] + sum(fq[k,t] for k in nodes[b].children) == fq[b,t])

            @constraint(m, v[b,t] <= nodes[b].v_max)
            @constraint(m, v[b,t] >= nodes[b].v_min)
        end
    end

    #for t in time_set
    #    @constraint(m, sum(rlp[b,t] for b in node_set) == sum(gp[b,t] + (wdp[b,t]-wdpc[b,t]) + (pvp[b,t]-pvpc[b,t]) - spc[b,t] + spd[b,t] for b in node_set))
    #    @constraint(m, sum(rlq[b,t] for b in node_set) == sum(gq[b,t] + wdq[b,t] + pvq[b,t] - sqc[b,t] + sqd[b,t] for b in node_set))
    #end

    # Voltage Equation of LinDistFlow for Multi-phase Systems

    #for b in setdiff(node_set, root_node)
    #    for t in time_set
    # All nodes without root
    #        b_ancestor = nodes[b].ancestor[1]
    #        @constraint(m, v[b,t] == v[b_ancestor,t] - 2*(lines_singlephase_to[b].r * fp[b,t] + lines_singlephase_to[b].x * fq[b,t]))
            #@constraint(m, lines_singlephase_to[b].s_max^2 >= fp[b]^2 + fq[b]^2)
    #    end
    #end

    for b in setdiff(node_set, root_node)
        for t in time_set
            b_ancestor = nodes[b].ancestor[1]
            bus_idx = nodes[b].bus
            bs = bus_idx
            samebus_node_idx = findall(bus_set .== bus_idx)
            samebusnodes = node_set[samebus_node_idx]
            samebusnodes_others = setdiff(samebusnodes, b)
            bb = samebusnodes_others
            if bb == []
                phase = []
                push!(phase, nodes[b].phase)
                @constraint(m, v[b,t] == v[b_ancestor,t] - 2*(lines_multiphase_to[bs].r[phase[1],phase[1]]*fp[b,t] + lines_multiphase_to[bs].x[phase[1],phase[1]]*fq[b,t]))
            else
                if length(bb) == 1
                    phases = []
                    push!(phases, nodes[b].phase)
                    for p in bb
                        push!(phases, nodes[p].phase)
                    end
                    if phases[1]==1 && phases[2]==2
                        #println("12")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                        - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                        - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+1,t]
                        + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+1,t])
                    elseif phases[1]==1 && phases[2]==3
                        #println("13")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                        - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                        + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+1,t]
                        - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+1,t])
                    elseif phases[1]==2 && phases[2]==1
                        #println("21")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                            - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                            + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-1,t]
                            - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b-1,t])
                    elseif phases[1]==2 && phases[2]==3
                        #println("23")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                        - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                        - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-1,t]
                        + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b-1,t])
                    elseif phases[1]==3 && phases[2]==1
                        #println("31")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                            - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                            - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-1,t]
                            + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b-1,t])

                    elseif phases[1]==3 && phases[2]==2
                        #println("32")
                        @constraint(m, v[b,t] == v[b_ancestor,t]
                            - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                            + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+1,t]
                            - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+1,t])
                    end
                else
                    phases = []
                    push!(phases, nodes[b].phase)
                    for p in bb
                        push!(phases, nodes[p].phase)
                    end
                    if phases[1]==1
                        #println("123")
                        if phases[2]==2
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+1,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+1,t]
                                    + (lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b+2,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]-lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b+2,t])
                        else
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b+1,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]+lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b+1,t]
                                    + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+2,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+2,t])
                        end
                    elseif phases[1]==2
                        #println("213")
                        if phases[2]==1
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-1,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[1]])*fq[b-1,t]
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b+1,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]+lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b+1,t])
                        else
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    + (lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b-1,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]-lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b-1,t]
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b+1,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b+1,t])
                        end
                    elseif phases[1]==3
                        #println("312")
                        if phases[2]==1
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-2,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]+lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b-2,t]
                                    + (lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b-1,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]-lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b-1,t])
                        else
                            @constraint(m, v[b,t] == v[b_ancestor,t]
                                    - 2*(lines_multiphase_to[bs].r[phases[1],phases[1]]*fp[b,t] + lines_multiphase_to[bs].x[phases[1],phases[1]]*fq[b,t])
                                    - (-lines_multiphase_to[bs].r[phases[1],phases[3]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[3]])*fp[b-2,t]
                                    + (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[3]]+lines_multiphase_to[bs].x[phases[1],phases[3]])*fq[b-2,t]
                                    + (lines_multiphase_to[bs].r[phases[1],phases[2]]+sqrt(3)*lines_multiphase_to[bs].x[phases[1],phases[2]])*fp[b-1,t]
                                    - (sqrt(3)*lines_multiphase_to[bs].r[phases[1],phases[2]]-lines_multiphase_to[bs].x[phases[1],phases[2]])*fq[b-1,t])
                        end
                    end
                end
            end
        end
    end

    #Voltage and flow constraints for root node(s)
    for b in root_node
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

    bus_results_P = DataFrame(bus = Any[], time = Any[], rl_P = Any[], g_P = Any[],
    wind_P = Any[], wind_Pcut = Any[], pv_P = Any[], pv_Pcut = Any[], s_Pch = Any[], s_Pdch = Any[], s_SOC = Any[], mu_PP = Any[])

    #bus_results_Q = DataFrame(bus = Any[], time = Any[], rl_Q = Any[], g_Q = Any[], s_Qch = Any[], s_Qdch = Any[], mu_QQ = Any[])
    bus_results_Q = DataFrame(bus = Any[], time = Any[], rl_Q = Any[], g_Q = Any[], wind_Q = Any[], pv_Q = Any[], s_Qch = Any[], s_Qdch = Any[], mu_QQ = Any[])

    bus_results_V = DataFrame(bus = Any[], time = Any[], voltage = Any[])

    line_results = DataFrame(line = Any[], from = Any[], to = Any[], time = Any[], f_P = Any[], f_Q = Any[], a_squared = Any[])

    for b in node_set
        for t in time_set
            resP = [b, t, getvalue(rlp[b,t])*Sbase/1000, getvalue(gp[b,t])*Sbase/1000, getvalue(wdp[b,t])*Sbase/1000,
            getvalue(wdpc[b,t])*Sbase/1000, getvalue(pvp[b,t])*Sbase/1000, getvalue(pvpc[b,t])*Sbase/1000, getvalue(spc[b,t])*Sbase/1000,
            getvalue(spd[b,t])*Sbase/1000, getvalue(ssoc[b,t])*100, getvalue(mup[b,t])]
            push!(bus_results_P, resP)

            resQ = [b, t, getvalue(rlq[b,t])*Sbase/1000, getvalue(gq[b,t])*Sbase/1000, getvalue(wdq[b,t])*Sbase/1000, getvalue(pvq[b,t])*Sbase/1000,
            getvalue(sqc[b,t])*Sbase/1000, getvalue(sqd[b,t])*Sbase/1000, getvalue(muq[b,t])]
            push!(bus_results_Q, resQ)
            resV = [b, t, sqrt(getvalue(v[b,t]))]
            push!(bus_results_V, resV)
        end
    end

    for b in node_set
        for t in time_set
            # Calc square current on line
            a_res = 0
            fp_res = getvalue(fp[b,t])*Sbase/1000
            fq_res = getvalue(fq[b,t])*Sbase/1000
            #if b != root_node
            if b in setdiff(node_set, root_node)
                v_res = getvalue(v[nodes[b].ancestor[1],t])
                a_res = (fp_res^2 + fq_res^2)/v_res
                lres = [lines_singlephase_to[b].index, nodes[b].ancestor[1], b, t, fp_res, fq_res, a_res]
                push!(line_results, lres)
            end
        end
    end

    sort!(bus_results_P)
    sort!(bus_results_Q)
    sort!(bus_results_V)
    sort!(line_results)

    gen_node = sort!(gen_node)
    wind_node = sort!(wind_node)
    pv_node = sort!(pv_node)
    stor_node = sort!(stor_node)

    P_restored = zeros(num_nodes, num_time_steps)
    mu_P = zeros(num_nodes, num_time_steps)
    for b in node_set
        P_restored[b,:] = bus_results_P[!,"rl_P"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
        #ArgumentError: syntax df[column] is not supported use df[!, column]
        mu_P[b,:] = bus_results_P[!,"mu_PP"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
    end

    Q_restored = zeros(num_nodes, num_time_steps)
    mu_Q = zeros(num_nodes, num_time_steps)
    for b in node_set
        Q_restored[b,:] = bus_results_Q[!,"rl_Q"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
        mu_Q[b,:] = bus_results_Q[!,"mu_QQ"][(b-1)*num_time_steps+1:(b-1)*num_time_steps+num_time_steps]
    end

    #Pmt_1p = zeros(length(gen_node), num_time_steps)
    #Qmt_1p = zeros(length(gen_node), num_time_steps)
    #m = 1
    #for b in gen_node
    #    Pmt_1p[m,:] = bus_results_P["g_P"][(gen_node[m]-1)*num_time_steps+1:(gen_node[m]-1)*num_time_steps+num_time_steps]
    #    Qmt_1p[m,:] = bus_results_Q["g_Q"][(gen_node[m]-1)*num_time_steps+1:(gen_node[m]-1)*num_time_steps+num_time_steps]
    #    m +=1
    #end
    #Pmt = zeros(1, num_time_steps)
    #Qmt = zeros(1, num_time_steps)
    #for mdx in 1:num_time_steps
    #    Pmt[mdx] = sum(Pmt_1p[:,mdx])
    #    Qmt[mdx] = sum(Qmt_1p[:,mdx])
    #end

    Pmt_ph1 = bus_results_P[!,"g_P"][(gen_node[1]-1)*num_time_steps+1:(gen_node[1]-1)*num_time_steps+num_time_steps]
    Pmt_ph2 = bus_results_P[!,"g_P"][(gen_node[2]-1)*num_time_steps+1:(gen_node[2]-1)*num_time_steps+num_time_steps]
    Pmt_ph3 = bus_results_P[!,"g_P"][(gen_node[3]-1)*num_time_steps+1:(gen_node[3]-1)*num_time_steps+num_time_steps]
    Pmt = Pmt_ph1 .+ Pmt_ph2 .+ Pmt_ph3
    Qmt_ph1 = bus_results_Q[!,"g_Q"][(gen_node[1]-1)*num_time_steps+1:(gen_node[1]-1)*num_time_steps+num_time_steps]
    Qmt_ph2 = bus_results_Q[!,"g_Q"][(gen_node[2]-1)*num_time_steps+1:(gen_node[2]-1)*num_time_steps+num_time_steps]
    Qmt_ph3 = bus_results_Q[!,"g_Q"][(gen_node[3]-1)*num_time_steps+1:(gen_node[3]-1)*num_time_steps+num_time_steps]
    Qmt = Qmt_ph1 .+ Qmt_ph2 .+ Qmt_ph3

    #Pwtb_1p = zeros(length(wind_node), num_time_steps)
    #Pwt_cut_1p = zeros(length(wind_node), num_time_steps)
    #w = 1
    #for b in wind_node
    #    Pwtb_1p[w,:] = bus_results_P["wind_P"][(wind_node[w]-1)*num_time_steps+1:(wind_node[w]-1)*num_time_steps+num_time_steps]
    #    Pwt_cut_1p[w,:] = bus_results_P["wind_Pcut"][(wind_node[w]-1)*num_time_steps+1:(wind_node[w]-1)*num_time_steps+num_time_steps]
    #    w +=1
    #end
    #Pwtb = zeros(1, num_time_steps)
    #Pwt_cut = zeros(1, num_time_steps)
    #for wdx in 1:num_time_steps
    #    Pwtb[wdx] = sum(Pwtb_1p[:,wdx])
    #    Pwt_cut[wdx] = sum(Pwt_cut_1p[:,wdx])
    #end

    Pwtb_ph1 = bus_results_P[!,"wind_P"][(wind_node[1]-1)*num_time_steps+1:(wind_node[1]-1)*num_time_steps+num_time_steps]
    Pwtb_ph2 = bus_results_P[!,"wind_P"][(wind_node[2]-1)*num_time_steps+1:(wind_node[2]-1)*num_time_steps+num_time_steps]
    Pwtb_ph3 = bus_results_P[!,"wind_P"][(wind_node[3]-1)*num_time_steps+1:(wind_node[3]-1)*num_time_steps+num_time_steps]
    Pwtb = Pwtb_ph1 .+ Pwtb_ph2 .+ Pwtb_ph3
    Pwt_cut_ph1 = bus_results_P[!,"wind_Pcut"][(wind_node[1]-1)*num_time_steps+1:(wind_node[1]-1)*num_time_steps+num_time_steps]
    Pwt_cut_ph2 = bus_results_P[!,"wind_Pcut"][(wind_node[2]-1)*num_time_steps+1:(wind_node[2]-1)*num_time_steps+num_time_steps]
    Pwt_cut_ph3 = bus_results_P[!,"wind_Pcut"][(wind_node[3]-1)*num_time_steps+1:(wind_node[3]-1)*num_time_steps+num_time_steps]
    Pwt_cut = Pwt_cut_ph1 .+ Pwt_cut_ph2 .+ Pwt_cut_ph3

    Qwt_inv_ph1 = bus_results_Q[!,"wind_Q"][(wind_node[1]-1)*num_time_steps+1:(wind_node[1]-1)*num_time_steps+num_time_steps]
    Qwt_inv_ph2 = bus_results_Q[!,"wind_Q"][(wind_node[2]-1)*num_time_steps+1:(wind_node[2]-1)*num_time_steps+num_time_steps]
    Qwt_inv_ph3 = bus_results_Q[!,"wind_Q"][(wind_node[3]-1)*num_time_steps+1:(wind_node[3]-1)*num_time_steps+num_time_steps]

    #Ppvs_1p = zeros(length(pv_node), num_time_steps)
    #Ppv_cut_1p = zeros(length(pv_node), num_time_steps)
    #p = 1
    #for b in pv_node
    #    Ppvs_1p[p,:] = bus_results_P["pv_P"][(pv_node[p]-1)*num_time_steps+1:(pv_node[p]-1)*num_time_steps+num_time_steps]
    #    Ppv_cut_1p[p,:] = bus_results_P["pv_Pcut"][(pv_node[p]-1)*num_time_steps+1:(pv_node[p]-1)*num_time_steps+num_time_steps]
    #    p +=1
    #end
    #Ppvs = zeros(1, num_time_steps)
    #Ppv_cut = zeros(1, num_time_steps)
    #for pdx in 1:num_time_steps
    #    Ppvs[pdx] = sum(Ppvs_1p[:,pdx])
    #    Ppv_cut[pdx] = sum(Ppv_cut_1p[:,pdx])
    #end

    Ppvs_ph1 = bus_results_P[!,"pv_P"][(pv_node[1]-1)*num_time_steps+1:(pv_node[1]-1)*num_time_steps+num_time_steps]
    Ppvs_ph2 = bus_results_P[!,"pv_P"][(pv_node[2]-1)*num_time_steps+1:(pv_node[2]-1)*num_time_steps+num_time_steps]
    Ppvs_ph3 = bus_results_P[!,"pv_P"][(pv_node[3]-1)*num_time_steps+1:(pv_node[3]-1)*num_time_steps+num_time_steps]
    Ppvs = Ppvs_ph1 .+ Ppvs_ph2 .+ Ppvs_ph3
    Ppv_cut_ph1 = bus_results_P[!,"pv_Pcut"][(pv_node[1]-1)*num_time_steps+1:(pv_node[1]-1)*num_time_steps+num_time_steps]
    Ppv_cut_ph2 = bus_results_P[!,"pv_Pcut"][(pv_node[2]-1)*num_time_steps+1:(pv_node[2]-1)*num_time_steps+num_time_steps]
    Ppv_cut_ph3 = bus_results_P[!,"pv_Pcut"][(pv_node[3]-1)*num_time_steps+1:(pv_node[3]-1)*num_time_steps+num_time_steps]
    Ppv_cut = Ppv_cut_ph1 .+ Ppv_cut_ph2 .+ Ppv_cut_ph3

    Qpv_inv_ph1 = bus_results_Q[!,"pv_Q"][(pv_node[1]-1)*num_time_steps+1:(pv_node[1]-1)*num_time_steps+num_time_steps]
    Qpv_inv_ph2 = bus_results_Q[!,"pv_Q"][(pv_node[2]-1)*num_time_steps+1:(pv_node[2]-1)*num_time_steps+num_time_steps]
    Qpv_inv_ph3 = bus_results_Q[!,"pv_Q"][(pv_node[3]-1)*num_time_steps+1:(pv_node[3]-1)*num_time_steps+num_time_steps]

    Pes_char1 = -1 .* bus_results_P[!,"s_Pch"][(stor_node[1]-1)*num_time_steps+1:(stor_node[1]-1)*num_time_steps+num_time_steps]
    Pes_dischar1 = bus_results_P[!,"s_Pdch"][(stor_node[1]-1)*num_time_steps+1:(stor_node[1]-1)*num_time_steps+num_time_steps]
    Pes1 = Pes_char1 + Pes_dischar1

    Pes_char2 = -1 .* bus_results_P[!,"s_Pch"][(stor_node[2]-1)*num_time_steps+1:(stor_node[2]-1)*num_time_steps+num_time_steps]
    Pes_dischar2 = bus_results_P[!,"s_Pdch"][(stor_node[2]-1)*num_time_steps+1:(stor_node[2]-1)*num_time_steps+num_time_steps]
    Pes2 = Pes_char2 + Pes_dischar2

    Pes_char3 = -1 .* bus_results_P[!,"s_Pch"][(stor_node[3]-1)*num_time_steps+1:(stor_node[3]-1)*num_time_steps+num_time_steps]
    Pes_dischar3 = bus_results_P[!,"s_Pdch"][(stor_node[3]-1)*num_time_steps+1:(stor_node[3]-1)*num_time_steps+num_time_steps]
    Pes3 = Pes_char3 + Pes_dischar3

    #Pes = Pes1 + Pes2 + Pes3
    #Pes = [Pes1 Pes2 Pes3]

    Qes_inv_cons1 = -1 .* bus_results_Q[!,"s_Qch"][(stor_node[1]-1)*num_time_steps+1:(stor_node[1]-1)*num_time_steps+num_time_steps]
    Qes_inv_inj1 = bus_results_Q[!,"s_Qdch"][(stor_node[1]-1)*num_time_steps+1:(stor_node[1]-1)*num_time_steps+num_time_steps]
    Qes1 = Qes_inv_cons1 + Qes_inv_inj1

    Qes_inv_cons2 = -1 .* bus_results_Q[!,"s_Qch"][(stor_node[2]-1)*num_time_steps+1:(stor_node[2]-1)*num_time_steps+num_time_steps]
    Qes_inv_inj2 = bus_results_Q[!,"s_Qdch"][(stor_node[2]-1)*num_time_steps+1:(stor_node[2]-1)*num_time_steps+num_time_steps]
    Qes2 = Qes_inv_cons2 + Qes_inv_inj2

    Qes_inv_cons3 = -1 .* bus_results_Q[!,"s_Qch"][(stor_node[3]-1)*num_time_steps+1:(stor_node[3]-1)*num_time_steps+num_time_steps]
    Qes_inv_inj3 = bus_results_Q[!,"s_Qdch"][(stor_node[3]-1)*num_time_steps+1:(stor_node[3]-1)*num_time_steps+num_time_steps]
    Qes3 = Qes_inv_cons3 + Qes_inv_inj3

    #Qes = Qes1 + Qes2 + Qes3
    #Qes = [Qes1 Qes2 Qes3]

    SOC_es1 = bus_results_P[!,"s_SOC"][(stor_node[1]-1)*num_time_steps+1:(stor_node[1]-1)*num_time_steps+num_time_steps]
    SOC_es2 = bus_results_P[!,"s_SOC"][(stor_node[2]-1)*num_time_steps+1:(stor_node[2]-1)*num_time_steps+num_time_steps]
    SOC_es3 = bus_results_P[!,"s_SOC"][(stor_node[3]-1)*num_time_steps+1:(stor_node[3]-1)*num_time_steps+num_time_steps]

    #SOC_es = (SOC_es1 + SOC_es2 + SOC_es3)/3
    #SOC_es = [SOC_es1 SOC_es2 SOC_es3]

    voltages = zeros(num_nodes, num_time_steps)
    for bdx in 1:num_nodes
        voltages[bdx,:] = bus_results_V[!,"voltage"][1+(bdx-1)*num_time_steps: bdx*num_time_steps]
    end

    fromnode = zeros(num_lines_singlephase)
    tonode = zeros(num_lines_singlephase)
    P_lineflow = zeros(num_lines_singlephase, num_time_steps)
    Q_lineflow = zeros(num_lines_singlephase, num_time_steps)
    for ldx in 1:num_lines_singlephase
        fromnode[ldx] = line_results[!,"from"][1+(ldx-1)*num_time_steps]
        tonode[ldx] = line_results[!,"to"][1+(ldx-1)*num_time_steps]
        P_lineflow[ldx,:] = line_results[!,"f_P"][1+(ldx-1)*num_time_steps:ldx*num_time_steps]
        Q_lineflow[ldx,:] = line_results[!,"f_Q"][1+(ldx-1)*num_time_steps:ldx*num_time_steps]
    end

    return objective_value, P_restored, Q_restored, Pmt, Qmt, Pwtb, Pwt_cut,
           Qwt_inv_ph1, Qwt_inv_ph2, Qwt_inv_ph3, Ppvs, Ppv_cut, Qpv_inv_ph1,
           Qpv_inv_ph2, Qpv_inv_ph3, Pes1, Pes2, Pes3, Qes1, Qes2, Qes3,
           SOC_es1, SOC_es2, SOC_es3, voltages, mu_P, mu_Q, fromnode, tonode, P_lineflow, Q_lineflow
end
