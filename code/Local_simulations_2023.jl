include("simulation_updated_2023.jl")
using Distributed 

world = Dict{Symbol, Any}(
    :nGens => 1,
    :realGen => 10,
    :worldSize => [3,3],
    :ratio => 0.5,
    :stab => 5,
    :fixed => [1,2],
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
