include("MongooseSimulation.jl")
using Distributed 
using Symbolics 
using DrWatson 
using StatsBase
using BSON

job_id = parse(Int, ARGS[1])
# job_id = 1

world = Dict{Symbol, Any}(
    :nGens => 500,
    :worldSize => [[3,3], [4,4], [5, 5]],
    :ratio => collect(0.05:0.05:0.95),
    # :ratio => 0.5,
    :stab => 2,
    :fixed => [[1,2]],
    :gain => 0.1,
    :loss => 0.1,
    :basem => [0.1, 0.01],
    :k => 0.1,
    :b => 1,
    :d => 0.1,
    :epsilon => 2,
    :multX => [0.01, 0.1],
    :multY => [0.01, 0.1],
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
out = produceOnceSim(out, false);

name = savename(cosm, "bson", accesses=cosm[:saveKeys])
location = joinpath("..", "data", name)
println("Saving, ", location)
wsave(location, out)
println(out[:err])
println("finished")

println("Running perturbations")
@time begin
out = re_evaluate_world_on_ratios(out, 0.1:0.1:0.9)
end
out[:ratio] = cosm[:ratio]
out[:peturbed] = true

name = savename(cosm, "bson", accesses=cosm[:saveKeys])
location = joinpath("..", "data", name)
println("Saving, ", location)
wsave(location, out)
println(out[:err])
println("finished")