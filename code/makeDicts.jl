using DrWatson 

world = Dict{Symbol, Any}(
    :force => [0.03, 0],
    :nGens => 1,
    :realGen =>  [@onlyif(:force != 0, 5000), @onlyif(:force==0, 10)],
    :q => 5,
    :n => 3,
    # :gain => [0.05, 0.1, 0.15, 0.2, 0.25],
    # :loss => [0.05, 0.1, 0.15, 0.2, 0.25],
    :stab => collect(5:1:15),
    :ratio => collect(0.1:0.1:0.9),
    # :ratio => vcat(0.1, collect(0.2:0.2:0.8), 0.9),
    :basem => 0.1,
    :k => 0.1,
    :b => 0.3,
    # :d => collect(0.1:0.1:0.9),    
    :d => [0.1, 0.5, 0.9, 
        @onlyif(:epsilon in (1, 5, 10), 
        [0.2, 0.3, 0.4, 0.6, 0.7, 0.8])...
    ],
    # :epsilon => collect(1:1:10),
    :epsilon => [1, 5, 10, 
        @onlyif(:d in (0.1, 0.5, 0.9), 
        [2, 3, 4, 6, 7, 8, 9])...
    ],
    :multX => 0.1,
    :multY => 0.1,
    :temp => Array{Any, 2}
)
world[:size] = world[:n]*world[:q]
worldSet = dict_list(world) 

for i in 1:length(worldSet)
    save(joinpath("tempDicts", string(i, ".bson")), worldSet[i])
end 

println(length(worldSet))