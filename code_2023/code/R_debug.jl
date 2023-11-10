include("MongooseSimulation.jl")

# for tRatio in 0.9:0.01:1

tRatio = 0.9

world = Dict{Symbol, Any}(
    :nGens => 500,
    :q => 3,
    :n => 3,
    :ratio => tRatio,
    :size => 9,
    # :ratio => 0.5,
    :stab => 1,
    :fixed => [[1,2]],
    :gain => 0.1,
    :loss => 0.1,
    :basem => 0.1,
    :k => 0.1,
    :b => 0.5,
    :d => 0.25,
    :epsilon => 0,
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
world[:gain] = world[:ratio]/world[:stab]
world[:loss] = (1-world[:ratio])/world[:stab]
world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))

@variables begin 
    F[1:world[:q], 1:world[:n]];
    W[1:world[:q], 1:world[:n]];
    R[1:world[:q], 1:world[:n]];
    M[1:world[:q], 1:world[:n]];
    Mf[1:world[:q], 1:world[:n]];
    Ml[1:world[:q], 1:world[:n]];
    P[1:world[:q], 1:world[:n]];
    Pf[1:world[:q], 1:world[:n]];
    Pl[1:world[:q], 1:world[:n]];
    C[1:world[:q], 1:world[:n]];
    Cf[1:world[:q], 1:world[:n]];
    Cl[1:world[:q], 1:world[:n]];
    Tr[1:world[:q], 1:world[:q]];
    X[1:world[:q], 1:world[:n]];
    Y[1:world[:q], 1:world[:n]];
    Xf[1:world[:q], 1:world[:n]];
    Yf[1:world[:q], 1:world[:n]];
    Xl[1:world[:q], 1:world[:n]];
    Yl[1:world[:q], 1:world[:n]];
end
# if world[:randStart]
#     world[:tX] = Symbolics.value.(rand(world[:q], world[:n]))
#     world[:tY] = Symbolics.value.(rand(world[:q], world[:n]))
# end
# world[:tX] = Symbolics.value.(rand(world[:q], world[:n])).*0.001
# world[:tX][:, 1] .= 0.0
# world[:tY] = Symbolics.value.(rand(world[:q], world[:n])).*0.001
# world[:tY][:, 1] .= 0.0

world[:update_x] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
world[:update_y] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
world[:cache_x] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
world[:cache_y] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))

u = ones(world[:q], world[:n])
u[:, 1] .= 0.0
world[:tW] = u
world[:tF] = reshape(
    repeat([1/world[:size]], world[:size]), (world[:q], world[:n])
)
world[:tR] = reshape(
    repeat(
        hcat([0 0], [1/p for p in 2:(world[:n]-1)]'), 
        world[:q]
    ), 
    (world[:q], world[:n])
)

modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC = makeModelExpr(world)

world = updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)

# create F system
fSys = makeFsys(F, M, P, C, d, world)
funF = build_function(
    fSys, F, M, P, C, d, epsilon;
    expression=Val{false}
    );
callFun = eval(funF[2]);

rSys = Symbolics.scalarize(makeRsys(R, M, P, C, d, world));
# R2 = collect(Symbolics.value.(R));
# [R2[p, 1] = 0.0 for p in 1:world[:q]];
# [R2[p, 2] = 0.0 for p in 1:world[:q]];
rSys = substitute(rSys, Dict(vcat([R[p, 1]=> 0 for p in 1:world[:q]], [R[p, 2]=> 0 for p in 1:world[:q]])))

funR = build_function(
    rSys, R, F, M, P, C, d, epsilon;
    expression=Val{false}
    );
callFunR = eval(funR[2]);

# create W system 
wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world);

# wSys[world[:fixed]...] = 1-W[world[:fixed]...]
# numInds = reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))'
# wSys[world[:fixed][1], world[:fixed][2]] = Symbolics.scalarize(1 - sum((W.*F.*numInds)./sum(F.*numInds)))

for q in 1:world[:q]
    wSys[q, 1] = W[q,1]
end
wSysSelec = Symbolics.scalarize(Wn./Wd);
funW = build_function(
    wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
    expression=Val{false}
    );
callFunW = eval(funW[2]);

grads, directSel, indirectSel = makeSelGrads(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, world);

fFun(Fx, x) = callFun(
    reshape(Fx, (world[:q], world[:n])), 
    reshape(x, (world[:q], world[:n])),  
    world[:M], 
    world[:P], 
    world[:C], 
    world[:d], 
    world[:epsilon]
    )

rFun(Fx, x) = callFunR(
    reshape(Fx, (world[:q], world[:n])), 
    reshape(x, (world[:q], world[:n])), 
    solF.zero, 
    world[:M], 
    world[:P], 
    world[:C], 
    world[:d], 
    world[:epsilon]
    )

wFun(Fx, x) = callFunW(
    reshape(Fx, (world[:q], world[:n])), 
    reshape(x, (world[:q], world[:n])), 
    solF.zero, 
    world[:Mf], 
    world[:Ml], 
    world[:Pf], 
    world[:Pl],
    world[:P],
    world[:C],
    world[:Cl],
    world[:Cf],
    world[:d], 
    world[:epsilon]
    )

# solF = zeros(world[:q], world[:n])
# try
#     try 
#         solF = genSolF2(fFun, world)
#     catch
#         solF = genSolF(fFun, world)
#     finally 
#         # println("F sol failed for:\n", 
#         #         savename(world))
#         # display(world[:tF])
#         # display(world[:tW])
        
#     end
# catch
#     # println("Defaulting to previous F")
#     # display(world[:tF])
#     solF = world[:tF]
# end


solF = genSolF2(fFun, world)
display(solF.zero)

solR = genSolR(rFun, world)
display(solR.zero)

outOFPlaceR = eval(funR[1])

function optimR(x)
    rp = zeros(world[:q], world[:n])
    rp[:, 3:end] .= reshape(x, (world[:q], world[:n]-2))
    return sum(abs.(callFunR(
        solF.zero, 
        rp, 
        world[:M], 
        world[:P], 
        world[:C], 
        world[:d], 
        world[:epsilon]
        )))
end

res = bboptimize(optimR, [0.5, 0.5, 0.5]; SearchRange = (0.0, 1.0), NumDimensions = 3, Method = :separable_nes)

outOFPlaceF = eval(funF[1])

function optimF(x)
    return sum(abs.(callFun(
        reshape(x, (world[:q], world[:n])), 
        world[:M], 
        world[:P], 
        world[:C], 
        world[:d], 
        world[:epsilon]
        )))
end
# end