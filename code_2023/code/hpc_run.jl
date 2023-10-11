using MongooseSimulation
using Distributed 
using Symbolics 
using DrWatson 
using StatsBase

# using BenchmarkTools

job_id = parse(Int, ARGS[1])

world = Dict{Symbol, Any}(
    :nGens => 250,
    :worldSize => [[5, 5]],
    :ratio => collect(0.05:0.05:0.95),
    # :ratio => 0.5,
    :stab => collect(2:0.25:3),
    :fixed => [[1,2]],
    :gain => 0.1,
    :loss => 0.1,
    :basem => 0.1,
    :k => 0.1,
    :b => 0.5,
    :d => 0.25,
    :epsilon => collect(2:0.25:3),
    :multX => 0.1,
    :multY => 0.1,
    :shape_X_cost => 2, 
    :shape_Y_cost => 2,
    :grad_rate => 0.1,
    :learning_rate => 0.01,
    :decay => 0.995,
    :verbose => true,
    :randStart => false,
    :SolFFails => 0,
    :SolWFails => 0,
    :SolRFails => 0,
    :saveKeys => [[
        :multX,
        :multY,
        :shape_X_cost,
        :shape_Y_cost,
        :b,
        :k,
        :d,
        :q,
        :n,
        :decay, 
        :epsilon, 
        :learning_rate, 
        :ratio, 
        :stab,]]
)

worldSet = dict_list(world);

cosm = worldSet[job_id];
# for cosm in worldSet
cosm[:gain] = cosm[:ratio]/cosm[:stab]
cosm[:loss] = (1-cosm[:ratio])/cosm[:stab]
cosm[:fix] = string(cosm[:fixed][1], cosm[:fixed][2])
cosm[:q] = cosm[:worldSize][1]
cosm[:n] = cosm[:worldSize][2] 
cosm[:size] = cosm[:n]*cosm[:q] 
cosm[:gradX] = zeros(cosm[:n], cosm[:q])
cosm[:gradY] = zeros(cosm[:n], cosm[:q])
cosm[:err] = 5
cosm[:err_list] = zeros(world[:nGens])
# end

@variables begin 
    F[1:cosm[:q], 1:cosm[:n]];
    W[1:cosm[:q], 1:cosm[:n]];
    R[1:cosm[:q], 1:cosm[:n]];
    M[1:cosm[:q], 1:cosm[:n]];
    Mf[1:cosm[:q], 1:cosm[:n]];
    Ml[1:cosm[:q], 1:cosm[:n]];
    P[1:cosm[:q], 1:cosm[:n]];
    Pf[1:cosm[:q], 1:cosm[:n]];
    Pl[1:cosm[:q], 1:cosm[:n]];
    C[1:cosm[:q], 1:cosm[:n]];
    Cf[1:cosm[:q], 1:cosm[:n]];
    Cl[1:cosm[:q], 1:cosm[:n]];
    Tr[1:cosm[:q], 1:cosm[:q]];
    X[1:cosm[:q], 1:cosm[:n]];
    Y[1:cosm[:q], 1:cosm[:n]];
    Xf[1:cosm[:q], 1:cosm[:n]];
    Yf[1:cosm[:q], 1:cosm[:n]];
    Xl[1:cosm[:q], 1:cosm[:n]];
    Yl[1:cosm[:q], 1:cosm[:n]];
end

println(string("Number of cosmologies: ", length(worldSet)))

println(cosm[:ratio], " --- ", cosm[:stab], " --- ", cosm[:epsilon])
# make temp dict and run function with 2 gens to compile runSim
compileCosm = copy(cosm)
compileCosm[:nGens] = 2
compileCosm[:verbose] = true
runSim(compileCosm)
# run actual sim with full nGens
out = copy(cosm);
out = produceOnceSim(out, true);
println(out[:err])

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

#  0.96 is expected at itr 40 for 6x6 should take ~30mins actually too 148s 
#  