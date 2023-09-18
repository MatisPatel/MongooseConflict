include("simulation_updated_2023.jl")
using Distributed 

world = Dict{Symbol, Any}(
    :nGens => 10,
    :worldSize => [[3,3]],
    :ratio => 0.25,
    :stab => 5,
    :fixed => [[1,2]],
    :q => 3,
    :n => 3,
    :gain => 0.1,
    :loss => 0.1,
    :basem => 0.1,
    :k => 0.1,
    :b => 0.3,
    :d => 0.1,
    :epsilon => collect(1:10),
    :multX => 0.1,
    :multY => 0.1,
    :force => 0.2
)
world[:q] = world[:worldSize][1][1] 
world[:n] = world[:worldSize][1][2] 
world[:size] = world[:n]*world[:q] 

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

@benchmark produceSim(worldSet[1], true)