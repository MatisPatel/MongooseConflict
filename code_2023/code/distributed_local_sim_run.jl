include("./simulation_updated_2023.jl")
using Distributed 
using BenchmarkTools
world = Dict{Symbol, Any}(
    :nGens => 150,
    :worldSize => [[5, 5]],
    :ratio => collect(0.1:0.1:0.9),
    # :ratio => 0.5,
    :stab => 2,
    :fixed => [[1,2]],
    :gain => 0.1,
    :loss => 0.1,
    :basem => 0.1,
    :k => [0.25],
    :b => [0.5],
    :d => [0.1],
    :epsilon => 2,
    :multX => 0.1,
    :multY => 0.1,
    :grad_rate => 0.1,
    :learning_rate => 0.01,
    :decay => 0.995,
    :verbose => false,
    :randStart => false,
    :SolFFails => 0,
    :SolWFails => 0,
    :SolRFails => 0,
)
world[:q] = world[:worldSize][1][1] 
world[:n] = world[:worldSize][1][2] 
world[:size] = world[:n]*world[:q] 
world[:gradX] = zeros(world[:n], world[:q])
world[:gradY] = zeros(world[:n], world[:q])
# world[:step_size_x] = zeros(world[:n], world[:q])
# world[:step_size_y] = zeros(world[:n], world[:q])
world[:err] = 5
world[:err_list] = [zeros(world[:nGens])]
# world[:step_size_x][:,2:end] .= world[:step_size_min]
# world[:step_size_y][:,2:end] .= world[:step_size_min]

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

println(string("Number of cosmologies: ", length(worldSet)))

Threads.@threads for cosm in worldSet
    println(cosm[:ratio], " --- ", cosm[:stab], " --- ", cosm[:epsilon])
    out = copy(cosm)
    out = produceOnceSim(out, true)
    println(out[:err])
end

# @time begin
# Threads.@threads for cosm in worldSet
#     println(Threads.threadid(), " --- ", cosm[:ratio], " --- ", cosm[:stab], " --- ", cosm[:epsilon])
#     out = copy(cosm)
#     for attempt in 1:5
#         out[:learning_rate] = cosm[:learning_rate]/(attempt)
#         println(out[:learning_rate])
#         out = produceSim(out, true)
#         println(out[:err])
#         if out[:err] < 1E-6
#             break
#         end
#     end
# end
# end

# decay_list = [0.995]
# lr_list = [0.01, 0.05, 0.005, 0.001]
# tt = length(decay_list)*length(lr_list)

# out_list = []

# for i in 1:length(decay_list)
#     for j in 1:length(lr_list)
#         worldSet[1][:decay] = decay_list[i]
#         worldSet[1][:learning_rate] = lr_list[j]
#         println("Decay: ", worldSet[1][:decay], "\nLR: ", worldSet[1][:learning_rate])
#         out = produceSim(worldSet[1], true)
#         push!(out_list, out)
#     end
# end

# using Plots
# plt = plot();
# for w in out_list
#     err = clamp.(w[:err_list][2:end], 0, 100)
#     plot!(plt, log1p.(err), label = string(w[:decay]," --- ", w[:learning_rate]));
# end

# save("test.png", plt)