# Data handler - Multi-phase Distribution Grid

using DataFrames
using CSV
using LinearAlgebra

# TYPE DEFINITIONS
mutable struct Generator
   index::Any
   bus::Int
   phase::Int
   node_idx::Int
   g_P_max::Float64
   g_Q_max::Float64
   cost::Float64
   function Generator(index, bus, phase, node_idx, g_P_max, g_S_max, cost)
      g = new()
      g.index  = index
      g.bus = bus
      g.phase = phase
      g.node_idx = node_idx
      g.g_P_max = g_P_max
      g.g_Q_max = sqrt(g_S_max^2 - g_P_max^2)
      g.cost = cost
      return g
   end
end

mutable struct Wind
   index::Any
   bus::Int
   phase::Int
   node_idx::Int
   w_P_max::Float64
   w_Q_max::Float64
   w_S_max::Float64
   function Wind(index, bus, phase, node_idx, w_P_max, w_S_max)
      w = new()
      w.index  = index
      w.bus = bus
      w.phase = phase
      w.node_idx = node_idx
      w.w_P_max = w_P_max
      w.w_S_max = w_S_max
      w.w_Q_max = sqrt(w_S_max^2 - w_P_max^2)
      return w
   end
end

mutable struct PV
   index::Any
   bus::Int
   phase::Int
   node_idx::Int
   p_P_max::Float64
   p_Q_max::Float64
   p_S_max::Float64
   function PV(index, bus, phase, node_idx, p_P_max, p_S_max)
      p = new()
      p.index  = index
      p.bus = bus
      p.phase = phase
      p.node_idx = node_idx
      p.p_P_max = p_P_max
      p.p_S_max = p_S_max
      p.p_Q_max = sqrt(p_S_max^2 - p_P_max^2)
      return p
   end
end

mutable struct Storage
   index::Any
   bus::Int
   phase::Int
   node_idx::Int
   s_P_max::Float64
   s_Q_max::Float64
   s_SOC_max::Float64
   s_SOC_min::Float64
   s_eff_char::Float64
   s_eff_dischar::Float64
   s_cap::Float64
   function Storage(index, bus, phase, node_idx, s_P_max, s_S_max, s_SOC_max, s_SOC_min, s_eff_char, s_eff_dischar, s_cap)
      s = new()
      s.index  = index
      s.bus = bus
      s.phase = phase
      s.node_idx = node_idx
      s.s_P_max = s_P_max
      s.s_Q_max = sqrt(s_S_max^2 - s_P_max^2)
      s.s_SOC_max = s_SOC_max
      s.s_SOC_min = s_SOC_min
      s.s_eff_char = s_eff_char
      s.s_eff_dischar = s_eff_dischar
      s.s_cap = s_cap
      return s
   end
end

mutable struct Node
   bus::Int
   phase::Int
   index::Any
   is_root::Bool
   d_P::Float64
   d_Q::Float64
   cosphi::Float64
   tanphi::Float64
   v_max::Float64
   v_min::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   generator::Generator
   wind::Wind
   pv::PV
   storage::Storage
   function Node(bus, phase, index, d_P, d_Q, v_max, v_min)
      b = new()
      b.bus = bus
      b.phase = phase
      b.index = index
      b.is_root = false
      b.d_P = d_P
      b.d_Q = d_Q
      b.v_max = v_max
      b.v_min = v_min
      b.children = Int[]
      b.ancestor = Int[]
      cosphi = d_P/(sqrt(d_P^2 + d_Q^2))
      tanphi = tan(acos(b.cosphi))
      if isnan(cosphi)
        b.cosphi = 0
        b.tanphi = 0
      else
        b.cosphi = cosphi
        b.tanphi = tan(acos(cosphi))
      end
      return b
   end
end

mutable struct LineSinglephase
   index::Any
   to_node::Int # the "to" node
   from_node::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   b::Float64 # the susceptance value
   s_max::Float64 # the capacity of the line
   function LineSinglephase(index, to_node, from_node, r, x, s_max)
      l = new()
      l.index = index
      l.to_node = to_node
      l.from_node = from_node
      l.r = r
      l.x = x
      l.b = (x/(r^2 + x^2))
      l.s_max = s_max
      return l
   end
end

mutable struct LineMultiphase
   index::Any
   to_bus::Int # the "to" bus
   from_bus::Int # the "from" bus
   r::Array{Float64} # the resistance matrix
   x::Array{Float64} # the reactance matrix
   b::Array{Float64} # the susceptance matrix
   s_max::Float64 # the capacity of the line
   function LineMultiphase(index, to_bus, from_bus, r, x, s_max)
      l = new()
      l.index = index
      l.to_bus = to_bus
      l.from_bus = from_bus
      l.r = r
      l.x = x
      l.b = (x./(r.^2 + x.^2))
      l.s_max = s_max
      return l
   end
end

function load_case_data(;datafile = "")

# READ RAW DATA

@info "Reading Data"

if datafile == ""
  data_dir = DATA_DIR
else
  data_dir = datafile
end

data_path = dirname(pwd())

nodes_raw = CSV.read(data_path*"/dataset/network_data/$data_dir/nodes_multiphase.csv", DataFrame)
sum(nonunique(nodes_raw, :index)) != 0 ? warn("Ambiguous Node Indices") : nothing

lines_raw_singlephase = CSV.read(data_path*"/dataset/network_data/$data_dir/lines_singlephase.csv", DataFrame)
sum(nonunique(lines_raw_singlephase, :index)) != 0  ? warn("Ambiguous Line Indices") : nothing

lines_raw_multiphase = CSV.read(data_path*"/dataset/network_data/$data_dir/lines_multiphase.csv", DataFrame)
sum(nonunique(lines_raw_multiphase, :index)) != 0  ? warn("Ambiguous Line Indices") : nothing

generators_raw = CSV.read(data_path*"/dataset/network_data/$data_dir/generators_multiphase.csv", DataFrame)
sum(nonunique(generators_raw, :index)) != 0 ? warn("Ambiguous Generator Indices") : nothing

windturbines_raw = CSV.read(data_path*"/dataset/network_data/$data_dir/windturbines_multiphase.csv", DataFrame)
sum(nonunique(windturbines_raw, :index)) != 0 ? warn("Ambiguous Wind Turbine Indices") : nothing

pvs_raw = CSV.read(data_path*"/dataset/network_data/$data_dir/pvs_multiphase.csv", DataFrame)
sum(nonunique(pvs_raw, :index)) != 0 ? warn("Ambiguous PV Indices") : nothing

storages_raw = CSV.read(data_path*"/dataset/network_data/$data_dir/storages_multiphase.csv", DataFrame)
sum(nonunique(storages_raw, :index)) != 0 ? warn("Ambiguous Storage Indices") : nothing

# Base values

Zbase = 1
Vbase = 2400
Sbase = (Vbase^2)/Zbase
Cbase = 800

# PREPARING MODEL DATA

nodes = Dict()
for n in 1:nrow(nodes_raw)
   bus = nodes_raw[n, :bus]
   phase = nodes_raw[n, :phase]
   index = nodes_raw[n, :index]
   d_P = 1000*nodes_raw[n, :d_P]/Sbase
   d_Q = 1000*nodes_raw[n, :d_Q]/Sbase
   v_max = nodes_raw[n, :v_max]
   v_min = nodes_raw[n, :v_min]
   newn = Node(bus, phase, index, d_P, d_Q, v_max, v_min)
   nodes[newn.index] = newn
end

lines_singlephase = Dict()
for l in 1:nrow(lines_raw_singlephase)
    index = lines_raw_singlephase[l, :index]
    from_node = lines_raw_singlephase[l, :from_node]
    to_node = lines_raw_singlephase[l, :to_node]
    r = lines_raw_singlephase[l, :r]/Zbase
    x = lines_raw_singlephase[l, :x]/Zbase
    s_max = 1000*lines_raw_singlephase[l, :s_max]/Sbase
    newl = LineSinglephase(index, to_node, from_node, r, x, s_max)

    push!(nodes[newl.from_node].children, newl.to_node)
    push!(nodes[newl.to_node].ancestor, newl.from_node)

    lines_singlephase[newl.index] = newl
end

lines_multiphase = Dict()
for l in 1:nrow(lines_raw_multiphase)
    index = lines_raw_multiphase[l, :index]
    from_bus = lines_raw_multiphase[l, :from_bus]
    to_bus = lines_raw_multiphase[l, :to_bus]

    # Line configurations - from the IEEE 13-bus system network data

    Zs1 = [0.3465+1.0179im 0.1560+0.5017im 0.1580+0.4236im;
             0.1560+0.5017im 0.3375+1.0478im 0.1535+0.3849im;
             0.1580+0.4236im 0.1535+0.3849im 0.3414+1.0348im]

     Ys1 = 10^-6*[6.2998im -1.9958im -1.2595im;
                   -1.9958im 5.9597im -0.7417im;
                   -1.2595im -0.7417im 5.6386im]

     # Configuration 602

     Zs2 = [0.7526+1.1814im 0.1580+0.4236im 0.1560+0.5017im;
             0.1580+0.4236im 0.7475+1.1983im 0.1535+0.3849im;
             0.1560+0.5017im 0.1535+0.3849im 0.7436+1.2112im]

     Ys2 = 10^-6*[5.6990im -1.0817im -1.6905im;
                   -1.0817im 5.1795im -0.6588im;
                   -1.6905im -0.6588im 5.4246im]


     # Configuration 603

     Zs3 = [0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
             0.0000+0.0000im 1.3294+1.3471im 0.2066+0.4591im;
             0.0000+0.0000im 0.2066+0.4591im 1.3238+1.3569im]

     Ys3 = 10^-6*[0.0000im 0.0000im 0.0000im;
                   0.0000im 4.7097im -0.8999im;
                   0.0000im -0.8999im 4.6658im]

     # Configuration 604

     Zs4 = [1.3238+1.3569im 0.0000+0.0000im 0.2066+0.4591im;
             0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
             0.2066+0.4591im 0.0000+0.0000im 1.3294+1.3471im]

     Ys4 = 10^-6*[4.6658im 0.0000im -0.8999im;
                   0.0000im 0.0000im 0.0000im;
                   -0.8999im 0.0000im 4.7097im]

     # Configuration 605

     Zs5 = [0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
             0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
             0.0000+0.0000im 0.0000+0.0000im 1.3292+1.3475im]

     Ys5 = 10^-6*[0.0000im 0.0000im 0.0000im;
                   0.0000im 0.0000im 0.0000im;
                   0.0000im 0.0000im 4.5193im]

     # Configuration 606

     Zs6 = [0.7982+0.4463im 0.3192+0.0328im 0.2849-0.0143im;
             0.3192+0.0328im 0.7891+0.4041im 0.3192+0.0328im;
             0.2849-0.0143im 0.3192+0.0328im 0.7982+0.4463im]

     Ys6 = 10^-6*[96.8897im 0.0000im 0.0000im;
                   0.0000im 96.8897im 0.0000im;
                   0.0000im 0.0000im 96.8897im]

     # Configuration 607

     Zs7 = [1.3425+0.5124im 0.0000+0.0000im 0.0000+0.0000im;
             0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im;
             0.0000+0.0000im 0.0000+0.0000im 0.0000+0.0000im]

     Ys7 = 10^-6*[88.9912im 0.0000im 0.0000im;
                   0.0000im 0.0000im 0.0000im;
                   0.0000im 0.0000im 0.0000im]

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

    # Update the lines_multiphase dict. with the line matrices

    if from_bus == 1 && to_bus == 4   # Line 650-632
      z = Z12
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 4 && to_bus == 7   # Line 632-670
      z = Z23
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 4 && to_bus == 5   # Line 632-633
      z = Z26
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 4 && to_bus == 3   # Line 632-645
      z = Z28
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 7 && to_bus == 10   # Line 670-671
      z = Z34
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 10 && to_bus == 14   # Line 671-680
      z = Z45
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 10 && to_bus == 11   # Line 671-692
      z = Z410
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 10 && to_bus == 9   # Line 671-684
      z = Z412
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 5 && to_bus == 6   # Line 633-634
      z = Z67
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 3 && to_bus == 2   # Line 645-646
      z = Z89
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 11 && to_bus == 12   # Line 692-675
      z = Z1011
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 9 && to_bus == 8   # Line 684-611
      z = Z1213
      r = real(z)/Zbase
      x = imag(z)/Zbase
   elseif from_bus == 9 && to_bus == 13   # Line 684-652
      z = Z1214
      r = real(z)/Zbase
      x = imag(z)/Zbase
   end

    s_max = 1000*lines_raw_multiphase[l, :s_max]/Sbase
    newl = LineMultiphase(index, to_bus, from_bus, r, x, s_max)

    lines_multiphase[newl.index] = newl
end

# Check topology
r = 0
root_node = 0
for b in keys(nodes)
    l = length(nodes[b].ancestor)
    if l > 1
        warn("Network not Radial (Node $(nodes[b].index))")
        #println("Network not Radial")
    elseif l == 0
        nodes[b].is_root = true
        root_node = b
        r += 1
    end
end
if r == 0
    warn("No root detected")
    root_node = 0
elseif r > 1
    #warn("More than one root detected")
    println("More than one root nodes are detected, probably the network is multi-phase")
end

generators = Dict()
for g in 1:nrow(generators_raw)
    index = generators_raw[g, :index]
    bus = generators_raw[g, :bus]
    phase = generators_raw[g, :phase]
    node_idx = generators_raw[g, :node]
    g_P_max = 1000*generators_raw[g, :p_max]/Sbase
    g_S_max = 1000*generators_raw[g, :s_max]/Sbase
    cost = generators_raw[g, :cost]
    newg = Generator(index, bus, phase, node_idx, g_P_max, g_S_max, cost)

    nodes[newg.node_idx].generator = newg

    generators[newg.index] = newg
end

windturbines = Dict()
for w in 1:nrow(windturbines_raw)
    index = windturbines_raw[w, :index]
    bus = windturbines_raw[w, :bus]
    phase = windturbines_raw[w, :phase]
    node_idx = windturbines_raw[w, :node]
    w_P_max = 1000*windturbines_raw[w, :p_max]/Sbase
    w_S_max = 1000*windturbines_raw[w, :s_max]/Sbase
    neww = Wind(index, bus, phase, node_idx, w_P_max, w_S_max)

    nodes[neww.node_idx].wind = neww

    windturbines[neww.index] = neww
end

pvs = Dict()
for p in 1:nrow(pvs_raw)
    index = pvs_raw[p, :index]
    bus = pvs_raw[p, :bus]
    phase = pvs_raw[p, :phase]
    node_idx = pvs_raw[p, :node]
    p_P_max = 1000*pvs_raw[p, :p_max]/Sbase
    p_S_max = 1000*pvs_raw[p, :s_max]/Sbase
    newp = PV(index, bus, phase, node_idx, p_P_max, p_S_max)

    nodes[newp.node_idx].pv = newp

    pvs[newp.index] = newp
end

storages = Dict()
for s in 1:nrow(storages_raw)
    index = storages_raw[s, :index]
    bus = storages_raw[s, :bus]
    phase = storages_raw[s, :phase]
    node_idx = storages_raw[s, :node]
    s_P_max = 1000*storages_raw[s, :p_max]/Sbase
    s_S_max = 1000*storages_raw[s, :s_max]/Sbase
    s_SOC_max = storages_raw[s, :SOC_max]
    s_SOC_min = storages_raw[s, :SOC_min]
    s_eff_char = storages_raw[s, :eff_char]
    s_eff_dischar = storages_raw[s, :eff_dischar]
    s_cap = storages_raw[s, :capacity]/Cbase
    news = Storage(index, bus, phase, node_idx, s_P_max, s_S_max, s_SOC_max, s_SOC_min, s_eff_char, s_eff_dischar, s_cap)

    nodes[news.node_idx].storage = news

    storages[news.index] = news
end

#info("Done preparing Data")
@info "Done preparing Data"
return nodes, lines_singlephase, lines_multiphase, generators, windturbines, pvs, storages
end
