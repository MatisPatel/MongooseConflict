using DrWatson 

world = Dict{Symbol, Any}(
    :force => [0.03],
    :nGens => 1,
    :realGen => 10000,
    # :realGen =>  [@onlyif(:force != 0, 10000), @onlyif(:force==0, 10)],
    :worldSize => [(4, 5)],
    # :q => 5,
    # :n => 3,
    # :gain => [0.05, 0.1, 0.15, 0.2, 0.25],
    # :loss => [0.05, 0.1, 0.15, 0.2, 0.25],
    :stab => [10],
    # :stab => 10,
    # :ratio => collect(0.1:0.1:0.9),
    :ratio => collect(Set(vcat(0.1, collect(0.2:0.2:0.8), 0.9, collect(0.3:0.025:0.7)))),
    :basem => 0.1,
    :k => 0.15,
    :b => 0.3,
    # :d => 0.5,
    :d => collect(0.1:0.1:0.9),    
    # :d => [0.1, 0.5, 0.9, 
    #     @onlyif(:epsilon in (1, 5, 10), 
    #     [0.2, 0.3, 0.4, 0.6, 0.7, 0.8])...
    # ],
    :epsilon => [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16],
    # :epsilon => collect(1:1:10),
    # :epsilon => [1, 5, 10, 
    #     @onlyif(:d in (0.1, 0.5, 0.9), 
    #     [2, 3, 4, 6, 7, 8, 9])...
    # ],
    :multX => 0.1,
    :multY => 0.1,
    :fixed => [[1, 2]],
    :temp => Array{Any, 2}
)
dict_list(Dict(:a=>0:0.1:10, :b=>0:1:5))
println("Making dict list")
worldSet = dict_list(world) 
# println("removing past files")
# foreach(rm, joinpath.("tempDicts", filter(endswith(".bson"), readdir("tempDicts"))))
println("save new dicts")
println(length(worldSet))
for i in 1:length(worldSet)
    if i%100 == 0
        println(i)
    end
    save(joinpath("tempDicts", string(i, ".bson")), worldSet[i])
end 
