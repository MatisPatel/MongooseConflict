include("MongooseSimulation.jl")
using Symbolics 

files = readdir("../data")
outfiles = readdir("../data2")

for file in files 

if file in outfiles
    continue
end

println("Processing file: ", file)

world = load("../data/$file")

println("Defining Vars ")
println("q: ", world[:q]) 
println("n: ", world[:n])
Symbolics.@variables begin 
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

evolved_X = copy(world[:tX])
evolved_Y = copy(world[:tY])
evolved_R = copy(world[:tR])
evolved_W = copy(world[:tW])
evolved_ratio = copy(world[:ratio])

modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC = makeModelExpr(world)

world = updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)

fSys = makeFsys(F, M, P, C, d, world)
funF = build_function(
    fSys, F, M, P, C, d, epsilon;
    expression=Val{false}
    );
callFun = eval(funF[2]);

rSys = Symbolics.scalarize(makeRsys(R, M, P, C, d, world));
rSys = substitute(rSys, Dict(vcat([R[p, 1]=> 0 for p in 1:world[:q]], [R[p, 2]=> 0 for p in 1:world[:q]])))

funR = build_function(
    rSys, R, F, M, P, C, d, epsilon;
    expression=Val{false}
    );
callFunR = eval(funR[2]);

wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world);

for q in 1:world[:q]
    wSys[q, 1] = W[q,1]
end
wSysSelec = Symbolics.scalarize(Wn./Wd);
funW = build_function(
    wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
    expression=Val{false}
    );
callFunW = eval(funW[2]);

fFun(Fx, x) = callFun(
    reshape(Fx, (world[:q], world[:n])), 
    reshape(x, (world[:q], world[:n])),  
    world[:M], 
    world[:P], 
    world[:C], 
    world[:d], 
    world[:epsilon]
    )

for ratio in 0.1:0.1:0.9
    # update world [:ratio] and rebuild equations and resolve them for tF, tW and tR 
    world[:ratio] = ratio
    world[:gain] = world[:ratio]/world[:stab]
    world[:loss] = (1-world[:ratio])/world[:stab]
    fSys = makeFsys(F, M, P, C, d, world)

    fSys = makeFsys(F, M, P, C, d, world)
    funF = build_function(
        fSys, F, M, P, C, d, epsilon;
        expression=Val{false}
        );
    callFun = eval(funF[2]);

    rSys = Symbolics.scalarize(makeRsys(R, M, P, C, d, world));
    rSys = substitute(rSys, Dict(vcat([R[p, 1]=> 0 for p in 1:world[:q]], [R[p, 2]=> 0 for p in 1:world[:q]])))

    funR = build_function(
        rSys, R, F, M, P, C, d, epsilon;
        expression=Val{false}
        );
    callFunR = eval(funR[2]);

    wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world);

    for q in 1:world[:q]
        wSys[q, 1] = W[q,1]
    end
    wSysSelec = Symbolics.scalarize(Wn./Wd);
    funW = build_function(
        wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
        expression=Val{false}
        );
    callFunW = eval(funW[2]);


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

    solF = zeros(world[:q], world[:n])
    try 
        solF = genSolF2(fFun, world)
    catch
        try 
            solF = genSolF(fFun, world)  
        catch
            solF = zeros(world[:q], world[:n])
        end
    end

    solR = zeros(world[:q], world[:n])
    try 
        solR = genSolR2(rFun, world)
    catch
        try
            solR = genSolR(rFun, world)  
        catch
            solR = zeros(world[:q], world[:n])
        end
    end

    solW = zeros(world[:q], world[:n])
    try 
        solW = genSolW2(wFun, world)
    catch
        try
            solW = genSolW(wFun, world)  
        catch 
            solW = zeros(world[:q], world[:n])
        end
    end

    # println("ratio: ", ratio)
    # println("evolved F: ")
    # display(world[:tF]) 
    # println("displaced F: ")
    # display(solF.zero)
    try 
        world[Symbol("tF", "_", ratio)] = solF.zero
    catch 
        world[Symbol("tF", "_", ratio)] = solF
    end

    # println("evolved R: ") 
    # display(world[:tR])
    # println("displaced R: ")
    # display(solR.zero)

    try
        world[Symbol("tR", "_", ratio)] = solR.zero
    catch
        world[Symbol("tR", "_", ratio)] = solR
    end

    # println("evolved W: ")
    # display(world[:tW])
    # println("displaced W: ")
    # display(solW.zero)
    try 
        world[Symbol("tW", "_", ratio)] = solW.zero
    catch
        world[Symbol("tW", "_", ratio)] = solW
    end
end

# remove old data file and save new one
save("../data2/$file", world)

end