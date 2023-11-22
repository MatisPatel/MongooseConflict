include("MongooseSimulation.jl")
using Distributed
using Symbolics 

files = readdir("../data")
outfiles = readdir("../data2")

Threads.@threads for file in files

if file in outfiles
    println("File already processed, skipping")
    continue
end

println("Processing file: ", file)

world = load("../data/$file")

world = re_evaluate_world_on_ratios(world, 0.1:0.1:0.9)

# save new one
save("../data2/$file", world)
end