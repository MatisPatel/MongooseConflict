include("./simulation_updated_2023.jl")
using Distributed 
using BenchmarkTools
world = Dict{Symbol, Any}(
    :nGens => 3,
    :worldSize => [[3,4]],
    :ratio => 0.25,
    :stab => 5,
    :fixed => [[1,2]],
    :gain => 0.1,
    :loss => 0.1,
    :basem => 0.1,
    :k => 0.1,
    :b => 0.3,
    :d => 0.1,
    :epsilon => 1,
    :multX => 0.1,
    :multY => 0.1,
    :step_size_max => 0.1, 
    :step_size_min => 1E7,
    :inc_factor => 1.2,
    :dec_factor => 0.5,
    :verbose => true,
)
world[:q] = world[:worldSize][1][1] 
world[:n] = world[:worldSize][1][2] 
world[:size] = world[:n]*world[:q] 
world[:gradX] = zeros(world[:n], world[:q])
world[:gradY] = zeros(world[:n], world[:q])
world[:step_size_x] = zeros(world[:n], world[:q])
world[:step_size_y] = zeros(world[:n], world[:q])
world[:step_size_x][:,2:end] .= world[:step_size_min]
world[:step_size_y][:,2:end] .= world[:step_size_min]

@variables begin 
    F[1:world[:q], 1:world[:n]]
    W[1:world[:q], 1:world[:n]]
    R[1:world[:q], 1:world[:n]]
    M[1:world[:q], 1:world[:n]]
    Mf[1:world[:q], 1:world[:n]]
    Ml[1:world[:q], 1:world[:n]]
    P[1:world[:q], 1:world[:n]]
    Pf[1:world[:q], 1:world[:n]]
    Pl[1:world[:q], 1:world[:n]]
    C[1:world[:q], 1:world[:n]]
    Cf[1:world[:q], 1:world[:n]]
    Cl[1:world[:q], 1:world[:n]]
    Tr[1:world[:q], 1:world[:q]]
    X[1:world[:q], 1:world[:n]]
    Y[1:world[:q], 1:world[:n]]
    Xf[1:world[:q], 1:world[:n]]
    Yf[1:world[:q], 1:world[:n]]
    Xl[1:world[:q], 1:world[:n]]
    Yl[1:world[:q], 1:world[:n]]
end

worldSet = dict_list(world)
for cosm in worldSet
    cosm[:gain] = cosm[:ratio]/cosm[:stab]
    cosm[:loss] = (1-cosm[:ratio])/cosm[:stab]
    cosm[:fix] = string(cosm[:fixed][1], cosm[:fixed][2])
end

# @time begin
# Threads.@threads for cosm in worldSet
#     produceSim(cosm, true)
# end
# end

out = produceSim(worldSet[1], true)