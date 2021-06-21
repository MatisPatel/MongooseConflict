using ModelingToolkit 
using NLsolve
using Symbolics
using DrWatson

@variables g, l, d, epsilon, k, b, B, deltaX, deltaY, dummy, t

world = load(joinpath("tempDicts", string(ENV["SLURM_ARRAY_TASK_ID"], ".bson")))

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


# world = Dict{Symbol, Any}(
#     :nGens => 500,
#     :q => 5,
#     :n => 3,
#     :gain => 0.1,
#     :loss => 0.1,
#     :basem => 0.1,
#     :k => 0.1,
#     :b => 0.3,
#     :d => 0.1,
#     :epsilon => 5,
#     :multX => 0.1,
#     :multY => 0.1,
#     :force => 0.2
# )

function simpleTrans(world)
    Tr = makeFArray(world[:q], world[:q])
    for i in 1:world[:q]
        for j in 1:world[:q]
            if (i+1) == j
                Tr[i, j] = world[:gain] 
            elseif (i-1) == j 
                Tr[i, j] = world[:loss] 
            end
        end
    end
    return Symbolics.value.(Tr) 
end

function matSub(pairs...)
    pairList = []
    for pair in pairs
        push!(pairList, vec(reshape(pair, (length(pair), 1))))
    end
    return (collect(Iterators.flatten((pairList))),)
end

function makeFArray(world)
    Fp = Array{Num}(undef, world[:q], world[:n])
    for q in 1:world[:q]
        for n in 1:world[:n]
            Fp[q,n] = 0
        end 
    end 
    return Fp 
end 

function makeFArray(Q, N)
    Fp = Array{Num}(undef, Q, N)
    for q in 1:Q
        for n in 1:N
            Fp[q,n] = 0
        end 
    end 
    return Fp 
end 

function addTrans(Fp, F, world)
    for q in 1:world[:q]
        for n in 1:world[:n]
            for newQ in 1:world[:q]
                if q > newQ
                    Fp[q, n] -= world[:loss] * F[q, n]
                    Fp[newQ, n] += world[:loss] * F[q, n]
                elseif q < newQ 
                    Fp[q, n] -= world[:gain] * F[q, n]
                    Fp[newQ, n] += world[:gain] * F[q, n]
                end
            end
        end 
    end 
    return Fp 
end

function addTrans(Wn, Wd, W, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            for newQ in 1:world[:q]
                if q > newQ
                    Wn[q, n] += world[:loss] * W[newQ, n]
                    Wd[q, n] += world[:loss]
                elseif newQ > q 
                    Wn[q, n] += world[:gain] * W[newQ, n]
                    Wd[q, n] += world[:gain]
                end
            end
        end 
    end 
    return Wn, Wd 
end

function addTransR(Rp, R, F, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            for newQ in 1:world[:q]
                if q > newQ
                    Rp[q, n] -= world[:loss] * F[q, n] * R[q, n]
                    Rp[newQ, n] += world[:loss] * F[q, n] * R[q, n]
                elseif newQ > q 
                    Rp[q, n] -= world[:gain] * F[q, n] * R[q, n]
                    Rp[newQ, n] += world[:gain] * F[q, n] * R[q, n]
                end
            end
        end
    end
    return Rp
end

function addMortality(Fp, F, M, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            Fp[q, n] -= M[q, n]*(n-1)*F[q, n]
            Fp[q, n-1] += M[q, n]*(n-1)*F[q, n]
        end 
    end 
    return Fp 
end

function addMortality(Wn, Wd, W, Mf, Ml, world)
    # focal dies
    for q in 1:world[:q]
        for n in 1:world[:n]
            Wn[q, n] += 0.0
            Wd[q, n] += Mf[q, n]
        end 
    end
    # if nonFocal dies
    for q in 1:world[:q]
        for n in 3:world[:n] 
            Wn[q, n] += (n-2) * Ml[q, n] * W[q, n-1]
            Wd[q, n] += (n-2) * Ml[q, n]
        end 
    end 
    return Wn, Wd
end

function addMortalityR(Rp, R, F, M, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            Rp[q, n] -= M[q, n] * (n-1) * F[q, n] * R[q, n]
            Rp[q, n-1] += M[q, n] * (n-1) * F[q, n] * R[q, n]
        end
    end
    return Rp 
end    

function addLocBirth(Fp, F, P, d, world)
    for q in 1:world[:q]
        for n in 1:world[:n]-1
            Fp[q, n] -= (n-1)*P[q, n]*(1-d)*F[q, n]
            Fp[q, n+1] += (n-1)*P[q, n]*(1-d)*F[q, n]
        end
    end
    return Fp
end

function addLocBirth(Wn, Wd, W, Pf, Pl, d, world)
    for q in 1:world[:q]
        for n in 2:world[:n]-1
            Wn[q, n] += (1-d) * Pf[q, n] * 2 * W[q, n+1]
            Wd[q, n] += (1-d) * Pf[q, n]  
        end
    end
    for q in 1:world[:q]
        for n in 3:world[:n]-1
            Wn[q, n] += (1-d) * (n-2) * Pl[q, n] * W[q, n+1]
            Wd[q, n] += (1-d) * (n-2) * Pl[q, n]  
        end
    end
    return Wn, Wd
end

function addLocBirthR(Rp, R, F, P, d, world)
    for q in 1:world[:q]
        for n in 2:(world[:n]-1)
            newN = n
            oldN = n-1
            Rp[q, n] -= oldN * P[q, n] * (1-d) * F[q, n] * R[q, n]
            Rp[q, n+1] += oldN * P[q, n] * (1-d) * F[q, n] *
            (((newN - 2)/newN)*R[q, n] + (2/newN)*
            ((1/oldN) + ((oldN-1)/oldN)*R[q,n]))
        end
    end
    return Rp
end

function addDistantBirth(Wn, Wd, W, F, Pf, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            for qNew in 1:world[:q]
                for nNew in 1:(world[:n]-1)
                    Wn[q, n] += d * Pf[q,n] * F[qNew, nNew] * (W[q,n] + W[qNew, nNew+1])
                    Wd[q, n] += d * Pf[q,n] * F[qNew, nNew]
                end
            end
        end
    end
    return Wn, Wd
end

function addImmigration(Fp, F, P, d, world)
    pbar = sum(F.*P * [n-1 for n in 1:world[:n]])
    for q in 1:world[:q]
        for n in 1:world[:n]-1
            Fp[q, n] -= d * pbar * F[q, n]
            Fp[q, n+1] += d * pbar * F[q, n]
        end
    end
    return Fp
end

function addImmigration(Wn, Wd, W, F, P, d, world)
    pbar = sum(F.*P * [n-1 for n in 1:world[:n]])
    for q in 1:world[:q]
        for n in 2:world[:n]-1
            Wn[q, n] += d * pbar * W[q, n+1]
            Wd[q, n] += d * pbar
        end
    end
    return Wn, Wd
end

function addImmigrationR(Rp, R, F, P, d, world)
    Pbar = sum(F.*P * [n-1 for n in 1:world[:n]])
    for q in 1:world[:q]
        for n in 2:world[:n]-1
            newN = n
            Rp[q, n] -= d * Pbar * F[q, n] * ((newN-2)/newN)*R[q, n]
            Rp[q, n+1] += d * Pbar * F[q, n] * ((newN-2)/newN)*R[q, n]            
        end
    end
    return Rp
end

function victory(a, b)
    return (a+0.000001)/(a+b+2*0.000001)
end

function addFights(Fp, F, C, epsilon, world)
    for q in 1:world[:q]
        for n in 1:world[:n]
            for qOpp in 1:world[:q]
                for nOpp in 1:world[:n]
                    # println("fight", F[q,n], " -> ", F[qOpp, nOpp])
                    # halved to not double count x->y and y->x
                    pEncounter  =  F[q, n] * F[qOpp, nOpp] * epsilon
                    pVictory    =  victory((n-1) * C[q, n], (nOpp-1) * C[qOpp, nOpp])
                    # print(pVictory)
                    fWin = pEncounter*pVictory
                    fLoss = pEncounter*(1-pVictory)
                    # fight logic 
                    # what is focal patch is minimal richness? then can only win losses dont change freq
                    # or if qOpp is maximum and q is less than it
                    if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
                        # println("only win")
                        Fp[q,n] -= fWin
                        Fp[q+1,n] += fWin
                    # what if focal patch is non-zero but other patch is 0. Then can only lose
                    elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
                        # println("only lose")
                        Fp[q,n] -= fLoss
                        Fp[q-1,n] += fLoss
                    # if equal richness and neither max not min then they can win or lose 
                    elseif (q == qOpp) & (q > 1) & (q < world[:q])
                        # println("both win or lose")
                        # q wins
                        Fp[q,n] -= fWin
                        Fp[q+1,n] += fWin 
                        # q loses
                        Fp[q,n] -= fLoss
                        Fp[q-1,n] += fLoss
                    # or if diffrent and not min or max
                    elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
                        # println("both win or lose")
                        # q wins
                        Fp[q,n] -= fWin
                        Fp[q+1,n] += fWin
                        # q loses
                        Fp[q,n] -= fLoss
                        Fp[q-1,n] += fLoss
                    # else
                        # println("NO MATCH") 
                    end
                end
            end
        end
    end                                    
    return Fp
end

function addFights(Wn, Wd, W, F, C, Cf, Cl, epsilon, world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            for qOpp in 1:world[:q]
                for nOpp in 1:world[:n]
                    pEncounter = F[q, n] * F[qOpp, nOpp] * epsilon
                    pVictory = victory(Cf[q, n] + (n-2) * Cl[q, n], (nOpp-1) * C[qOpp, nOpp])
                    fWin = pEncounter * pVictory
                    fLoss = pEncounter * (1 - pVictory)
                    if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
                        # println("only win")
                        Wn[q, n] += fWin * W[q + 1, n]
                        Wd[q, n] += fWin
                    # what if focal patch is non-zero but other patch is 0. Then can only lose
                    elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
                        # println("only lose")
                        Wn[q, n] += fLoss * W[q - 1, n]
                        Wd[q, n] += fLoss
                    # if equal richness and neither max not min then they can win or lose 
                    elseif (q == qOpp) & (q > 1) & (q < world[:q])
                        # println("both win or lose")
                        # q wins
                        Wn[q, n] += fWin * W[q + 1, n]
                        Wd[q, n] += fWin 
                        # q loses
                        Wn[q, n] += fLoss * W[q - 1, n]
                        Wd[q, n] += fLoss
                    # or if diffrent and not min or max
                    elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
                        # println("both win or lose")
                        Wn[q, n] += fWin * W[q + 1, n]
                        Wd[q, n] += fWin 
                        # q loses
                        Wn[q, n] += fLoss * W[q - 1, n]
                        Wd[q, n] += fLoss
                    # else
                        # println("NO MATCH") 
                    end
                end
            end
        end
    end
    return Wn, Wd
end

function addFightsR(Rp, R, F, C, epsilon, world)
    for q in 1:world[:q]
        for n in 1:world[:n]
            for qOpp in 1:world[:q]
                for nOpp in 1:world[:n]
                    # println("fight", F[q,n], " -> ", F[qOpp, nOpp])
                    # halved to not double count x->y and y->x
                    pEncounter  =  F[q, n] * F[qOpp, nOpp] * epsilon
                    pVictory    =  victory((n-1) * C[q, n], (nOpp-1) * C[qOpp, nOpp])
                    # print(pVictory)
                    fWin = pEncounter*pVictory
                    fLoss = pEncounter*(1-pVictory)
                    # fight logic 
                    # what is focal patch is minimal richness? then can only win losses dont change freq
                    # or if qOpp is maximum and q is less than it
                    if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
                        # println("only win")
                        Rp[q,n] -= fWin * R[q, n]
                        Rp[q+1,n] += fWin * R[q, n]
                    # what if focal patch is non-zero but other patch is 0. Then can only lose
                    elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
                        # println("only lose")
                        Rp[q,n] -= fLoss * R[q,n]
                        Rp[q-1,n] += fLoss * R[q,n]
                    # if equal richness and neither max not min then they can win or lose 
                    elseif (q == qOpp) & (q > 1) & (q < world[:q])
                        # println("both win or lose")
                        # q wins
                        Rp[q,n] -= fWin * R[q,n]
                        Rp[q+1,n] += fWin * R[q,n]
                        # q loses
                        Rp[q,n] -= fLoss * R[q,n]
                        Rp[q-1,n] += fLoss * R[q,n]
                    # or if diffrent and not min or max
                    elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
                        # println("both win or lose")
                        # q wins
                        Rp[q,n] -= fWin * R[q,n]
                        Rp[q+1,n] += fWin * R[q,n]
                        # q loses
                        Rp[q,n] -= fLoss * R[q,n]
                        Rp[q-1,n] += fLoss * R[q,n]
                    # else
                        # println("NO MATCH") 
                    end
                end
            end
        end
    end                                    
    return Rp
end

function makeFsys(F, M, P, C, d, world)
    Fp = makeFArray(world)
    Fp1 = addTrans(copy(Fp), F, world)
    Fp2 = addMortality(copy(Fp1), F, M, world)
    Fp3 = addLocBirth(copy(Fp2), F, P, d, world)
    Fp4 = addImmigration(copy(Fp3), F, P, d, world)
    Fp5 = addFights(copy(Fp4), F, C, epsilon, world)
    Fp5[1,2] = 1 - sum(F)
    return Fp5
end

function mortFun2(n, xf, xl, yf, yl, world)
    calc = B * 2.718^(-1 * (n-1)*xl) + world[:multX]*xf^2 + world[:multY]*yf^2
    return calc
end 

function mortFun(n, xf, xl, yf, yl, world)
    calc = B * 2.718^(-1 * (n)*((xl*(n-1) + xf)/n)) + world[:multX]*xf^2 + world[:multY]*yf^2
    return calc
end 

function makeFocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = mortFun(n, Xf[q, n], Xl[q, n], Yf[q, n], Yl[q, n], world)
        end
    end
    return model
end

function makeLocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = mortFun(n, Xl[q, n], Xf[q, n], Yl[q, n], Yf[q, n], world)
        end
    end
    return model
end

function prodFun(q, n, world)
    calc = world[:k] + ((q-1)*world[:b])*(1/(n-1))
    return calc
end

function makeFocalModelP(world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = prodFun(q, n, world)
        end
        model[q, 1] = 0.01
    end
    return model
end

function makeLocalModelP(world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 3:world[:n]
            model[q, n] = prodFun(q, n, world)
        end
        model[q, 1] = 0.01
    end
    return model
end

function makeFocalModelC(Yf, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = Yf[q, n]
        end
    end
    return model
end

function makeLocalModelC(Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 3:world[:n]
            model[q, n] = Yl[q, n]
        end
    end
    return model
end

function makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)
    Wn = makeFArray(world)
    Wd = makeFArray(world)
    Wn1, Wd1 = addMortality(copy(Wn), copy(Wd), W, Mf, Ml, world)
    Wn2, Wd2 = addTrans(copy(Wn1), copy(Wd1), W, world)
    Wn3, Wd3 = addLocBirth(copy(Wn2), copy(Wd2), W, Pf, Pl, d, world)
    Wn4, Wd4 = addImmigration(copy(Wn3), copy(Wd3), W, F, P, d, world)
    Wn5, Wd5 = addFights(copy(Wn4), copy(Wd4), W, F, C, Cf, Cl, epsilon, world)
    Wn6, Wd6 = addDistantBirth(copy(Wn5), copy(Wd5), W, F, Pf, world)
    wSys = Wn6./Wd6.-W
    # wSys[1,2] = 1-W[1,2]
    # for q in 1:world[:q]
    #     wSys[q, 1] = W[q,1]
    # end
    return wSys, Wn6, Wd6 
end

function makeRsys(R, M, P, C, d, world)
    Rp = makeFArray(world)
    Rp1 = addTransR(copy(Rp), R, F, world)
    Rp2 = addMortalityR(copy(Rp1), R, F, M, world)
    Rp3 = addLocBirthR(copy(Rp2), R, F, P, d, world)
    Rp4 = addImmigrationR(copy(Rp3), R, F, P, d, world)
    Rp5 = addFightsR(copy(Rp4), R, F, C, epsilon, world)
    return Rp5
end

atEqXY = matSub(Xf .=> X, Xl .=> X, Yf .=> Y, Yl .=> Y)
atEqMPC = matSub(Mf .=> M, Ml .=> M, Pf .=> P, Pl .=> P, Cf .=> C, Cl .=> C)

function makeModelExpr(world)
    modelMf = makeFocalModelM(Xf, Xl, Yf, Yl, world)
    modelMl = makeLocalModelM(Xf, Xl, Yf, Yl, world)

    modelPf = makeFocalModelP(world)
    modelPl = makeLocalModelP(world)

    modelCf = makeFocalModelC(Yf, world)
    modelCl = makeLocalModelC(Yl, world)

    atEqXY = matSub(Xf .=> X, Xl .=> X, Yf .=> Y, Yl .=> Y)
    modelM = substitute.(modelMf, atEqXY)
    modelP = substitute.(modelPf, atEqXY)
    modelC = substitute.(modelCf, atEqXY)

    modelMFun = eval(ModelingToolkit.build_function(
        modelM, B, X, Y
    )[1])
    modelPFun = eval(ModelingToolkit.build_function(
        modelP, b, k
    )[1])
    modelCFun = eval(ModelingToolkit.build_function(
        modelC, Y
    )[1])

    modelMlFun = eval(ModelingToolkit.build_function(
        modelM, B, Xf, Xl, Yl
    )[1])
    modelPlFun = eval(ModelingToolkit.build_function(
        modelP, b, k
    )[1])
    modelClFun = eval(ModelingToolkit.build_function(
        modelC, Yl
    )[1])
    return modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC
end

# function updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)
#     world[:M] = modelMFun(world[:basem], world[:tX], world[:tY])
#     world[:P] = modelPFun(world[:b], world[:k])
#     world[:C] = modelCFun(world[:tY])
#     world[:Mf] = world[:M]
#     world[:Pf] = world[:P]
#     world[:Cf] = world[:C]
#     world[:Ml] = modelMlFun(world[:basem], world[:tX], world[:tX], world[:tY])
#     world[:Pl] = modelPlFun(world[:b], world[:k])
#     world[:Cl] = modelClFun(world[:tY])
#     return world
# end

function updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)
    world[:M] = Symbolics.value.(substitute.(modelM, (vcat([B=>world[:basem]], matSub(X.=>world[:tX], Y.=>world[:tY])[1]),))) 
    world[:P] = Symbolics.value.(modelP)
    world[:C] = Symbolics.value.(substitute.(modelC, matSub(Y.=>world[:tY])))
    world[:Mf] = world[:M]
    world[:Pf] = world[:P]
    world[:Cf] = world[:C]
    world[:Ml] = Symbolics.value.(substitute.(modelMl, 
        (vcat([B=>world[:basem]], 
        matSub(
            Xf.=>world[:tX], 
            Xl.=>world[:tX], 
            Yl.=>world[:tY],
            )[1]),)
    ))
    world[:Pl] = Symbolics.value.(modelPl)
    world[:Cl] = Symbolics.value.(substitute.(modelCl, matSub(Yl.=>world[:tY])))
    return world
end

function evalDeriv(Fun::Array{Num}, var, sub)
    return substitute.(Symbolics.derivative.(Fun, var), (sub,))
end

function evalDeriv(Fun::Num, var, sub)
    return substitute.(Symbolics.derivative.(Fun, var), (sub,))[1]
end

function makeSelGrads(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, 
    modelPf, modelPl, modelCf, modelCl, world)
    atEqMPC = matSub(Mf .=> M, Ml .=> M, Pf .=> P, Pl .=> P, Cf .=> C, Cl .=> C)

    gradsMf = Symbolics.derivative.(wSysSelec, Mf)
    gradsPf = Symbolics.derivative.(wSysSelec, Pf)
    gradsCf = Symbolics.derivative.(wSysSelec, Cf)
    gradsMl = Symbolics.derivative.(wSysSelec, Ml)
    gradsPl = Symbolics.derivative.(wSysSelec, Pl)
    gradsCl = Symbolics.derivative.(wSysSelec, Cl)

    gradsF = [gradsMf, gradsPf, gradsCf]
    gradsFEq = [substitute.(x, atEqMPC) for x in gradsF]

    gradsL = [gradsMl, gradsPl, gradsCl]
    gradsLEq = [substitute.(x, atEqMPC) for x in gradsL]

    # args = vcat(matSub(W.=>world[:tW], F.=>world[:tF], Tr.=>tTr, M.=>tM, P.=>tP, C.=>tC)[1], [d=>world[:d], epsilon=>world[:epsilon]])

    # selGF = [substitute.(x, (args,)) for x in gradsFEq]
    # selGL = [substitute.(x, (args,)) for x in gradsLEq]
    selGF = gradsFEq
    selGL = gradsLEq

    # eqs = reshape(0 .~ Fp4, (6,1))
    # sts = reshape(F, (6,1))
    # pars = [Tr, M, P, d]
    atEqXY = matSub(Xf.=>X, Xl.=>X, Yf.=>Y, Yl.=>Y)

    wSub = substitute.(wSysSelec, matSub(Mf.=>modelMf, Pf.=>modelPf, Cf.=>modelCf, Ml.=>modelMl, Pl.=>modelPl, Cl.=>modelCl))
    wDiffXf = Symbolics.derivative.(wSub, Xf)
    wDiffXl = Symbolics.derivative.(wSub, Xl)
    wDiffYf = Symbolics.derivative.(wSub, Yf)
    wDiffYl = Symbolics.derivative.(wSub, Yl)

    nX = makeFArray(world)
    nY = copy(nX)
    nXDirect = copy(nX)
    nYDirect = copy(nX)
    nXInd = copy(nX) 
    nYInd = copy(nY)

    for q in 1:world[:q]
        for n in 1:world[:n]
            directX = 
            evalDeriv(modelMf[q, n], Xf[q, n], atEqXY[1])*selGF[1][q, n] +
            evalDeriv(modelMl[q, n], Xf[q, n], atEqXY[1])*selGL[1][q, n] +
            evalDeriv(modelPf[q, n], Xf[q, n], atEqXY[1])*selGF[2][q, n] +
            evalDeriv(modelPl[q, n], Xf[q, n], atEqXY[1])*selGL[2][q, n] +
            evalDeriv(modelCf[q, n], Xf[q, n], atEqXY[1])*selGF[3][q, n] +
            evalDeriv(modelCl[q, n], Xf[q, n], atEqXY[1])*selGL[3][q, n]

            indirectX = 
            evalDeriv(modelMf[q, n], Xl[q, n], atEqXY[1])*selGF[1][q, n] +
            evalDeriv(modelMl[q, n], Xl[q, n], atEqXY[1])*selGL[1][q, n] +
            evalDeriv(modelPf[q, n], Xl[q, n], atEqXY[1])*selGF[2][q, n] +
            evalDeriv(modelPl[q, n], Xl[q, n], atEqXY[1])*selGL[2][q, n] +
            evalDeriv(modelCf[q, n], Xl[q, n], atEqXY[1])*selGF[3][q, n] +
            evalDeriv(modelCl[q, n], Xl[q, n], atEqXY[1])*selGL[3][q, n]

            directY = 
            evalDeriv(modelMf[q, n], Yf[q, n], atEqXY[1])*selGF[1][q, n] +
            evalDeriv(modelMl[q, n], Yf[q, n], atEqXY[1])*selGL[1][q, n] +
            evalDeriv(modelPf[q, n], Yf[q, n], atEqXY[1])*selGF[2][q, n] +
            evalDeriv(modelPl[q, n], Yf[q, n], atEqXY[1])*selGL[2][q, n] +
            evalDeriv(modelCf[q, n], Yf[q, n], atEqXY[1])*selGF[3][q, n] +
            evalDeriv(modelCl[q, n], Yf[q, n], atEqXY[1])*selGL[3][q, n]

            indirectY = 
            evalDeriv(modelMf[q, n], Yl[q, n], atEqXY[1])*selGF[1][q, n] +
            evalDeriv(modelMl[q, n], Yl[q, n], atEqXY[1])*selGL[1][q, n] +
            evalDeriv(modelPf[q, n], Yl[q, n], atEqXY[1])*selGF[2][q, n] +
            evalDeriv(modelPl[q, n], Yl[q, n], atEqXY[1])*selGL[2][q, n] +
            evalDeriv(modelCf[q, n], Yl[q, n], atEqXY[1])*selGF[3][q, n] +
            evalDeriv(modelCl[q, n], Yl[q, n], atEqXY[1])*selGL[3][q, n]

            nX[q, n] = directX + R[q, n] * indirectX 
            nXDirect[q, n] = directX 
            nXInd[q, n] = indirectX 

            nY[q, n] = directY + R[q, n] * indirectY 
            nYDirect[q, n] = directY 
            nYInd[q, n] = indirectY 
        end
    end

    # grads = [nX, nY]
    # DirectSel = [nXDirect, nYDirect]
    # IndirectSel = [nXInd, nYInd]
    subs = matSub(
            M.=>modelM, 
            P.=>modelP, 
            C.=>modelC
        )

    subsX = matSub(Xf.=>X, Xl.=>X, Yf.=>Y, Yl.=>Y)

    nX2 = substitute.(substitute.(nX, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))
    nY2 = substitute.(substitute.(nY, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))
    nXDirect2 = substitute.(substitute.(nXDirect, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))
    nYDirect2 = substitute.(substitute.(nYDirect, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))
    nXInd2 = substitute.(substitute.(nXInd, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))
    nYInd2 = substitute.(substitute.(nYInd, subs), 
        ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],))

    grads = [nX2, nY2]
    directSel = [nXDirect2, nYDirect2]
    indirectSel = [nXInd2, nYInd2]

    return grads, directSel, indirectSel
end

function genSolF(fFun, world)
    return nlsolve(fFun, 
        reshape(
            repeat([1/world[:size]], world[:size]), (world[:q], world[:n])
        )
    )
end

function genSolF2(fFun, world)
    return nlsolve(fFun, world[:tF])
end

function genSolW(wFun, world)
    u = ones(world[:q], world[:n])
    u[:, 1] .= 0.0
    return nlsolve(wFun, u)
end

function genSolW2(wFun, world)
    # u = ones(world[:q], world[:n])
    # u[:, 1] .= 0.0
    return nlsolve(wFun, world[:tW])
end

function genSolR2(rFun, world)
    return nlsolve(
        rFun, world[:tR]
        )
end

function genSolR(rFun, world)
    return nlsolve(
        rFun,
            reshape(
                repeat(
                    vcat([0, 0], [1/p for p in 2:(world[:n]-1)]), 
                    world[:q]
                ), 
                (world[:q], world[:n])
            )
        )
end

# world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 0.0))
# world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 0.0))
# world[:tTr] = simpleTrans(world)

# modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC = makeModelExpr(world)

# world = updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)

# # create F system
# fSys = makeFsys(F, M, P, C, d, world)
# funF = ModelingToolkit.build_function(
#     fSys, F, M, P, C, d, epsilon;
#     expression=Val{false}
#     );
# callFun = eval(funF[2]);

# rSys = makeRsys(R, M, P, C, d, world)
# R2 = copy(R);
# [R2[p, 1] = 0 for p in 1:world[:q]];
# [R2[p, 2] = 0 for p in 1:world[:q]];
# rSys = substitute.(rSys, matSub(R.=>R2))
# for q in 1:world[:q]
#     for n in 1:2
#         rSys[q, n] = R[q, n]
#     end 
# end 
# funR = ModelingToolkit.build_function(
#     rSys, R, F, M, P, C, d, epsilon;
#     expression=Val{false}
#     );
# callFunR = eval(funR[2]);

# # create W system 
# wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)
# wSys[2,1] = 1-W[2,1]
# for q in 1:world[:q]
#     wSys[q, 1] = W[q,1]
# end
# wSysSelec = Wn./Wd
# funW = ModelingToolkit.build_function(
#     wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
#     expression=Val{false}
#     );
# callFunW = eval(funW[2]);

# grads, directSel, indirectSel = makeSelGrads(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, world);
function corrErr(g, t)
    s = sign.(g)
    mask = g[((t.==0) .& (s .> 0)) .| ((t.==1) .& (s .< 0)) .| ((t .!= 1) .& (t .!= 0))]
    if isempty(mask)
        return 0.0
    else
        return sum(abs.(mask))
    end
end

function step(world, callFun, callFunW, callFunR, 
    grads, 
    modelM, modelP, modelC, modelMl, modelPl, modelCl
    )
    world = updateMPC(copy(world), modelM, modelP, modelC, modelMl, modelPl, modelCl)

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
        world[:Cf],
        world[:Cl],
        world[:d], 
        world[:epsilon]
        )

    
    solF = genSolF2(fFun, world)
    world[:tF] = solF.zero

    solW = genSolW2(wFun, world)
    world[:tW] = solW.zero

    solR = genSolR2(rFun, world)
    world[:tR] = solR.zero

    gradX = substitute.(grads[1], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

    gradY = substitute.(grads[2], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

    world[:gradX] = Symbolics.value.(gradX)
    world[:gradY] = Symbolics.value.(gradY)

    world[:tX] = clamp.(
        world[:tX] .+ world[:force]*clamp.(world[:gradX], -1, 1),
        # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
        0,
        1
    )

    world[:tY] = clamp.(
        world[:tY] .+ world[:force]*clamp.(world[:gradY], -1, 1),
        # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
        0,
        1
    )

    return world
end

function runSim(world)
    world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 0.0))
    world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 0.0))
    world[:tTr] = simpleTrans(world)

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
    funF = ModelingToolkit.build_function(
        fSys, F, M, P, C, d, epsilon;
        expression=Val{false}
        );
    callFun = eval(funF[2]);

    rSys = makeRsys(R, M, P, C, d, world)
    R2 = copy(R);
    [R2[p, 1] = 0 for p in 1:world[:q]];
    [R2[p, 2] = 0 for p in 1:world[:q]];
    rSys = substitute.(rSys, matSub(R.=>R2))
    for q in 1:world[:q]
        for n in 1:2
            rSys[q, n] = R[q, n]
        end 
    end 
    funR = ModelingToolkit.build_function(
        rSys, R, F, M, P, C, d, epsilon;
        expression=Val{false}
        );
    callFunR = eval(funR[2]);

    # create W system 
    wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)
    wSys[world[:fixed]...] = 1-W[world[:fixed]...]
    for q in 1:world[:q]
        wSys[q, 1] = W[q,1]
    end
    wSysSelec = Wn./Wd
    funW = ModelingToolkit.build_function(
        wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
        expression=Val{false}
        );
    callFunW = eval(funW[2]);

    grads, directSel, indirectSel = makeSelGrads(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, world);

    for i in 1:world[:nGens]
        world = step(world, callFun, callFunW, callFunR,
            grads,
            modelM, modelP, modelC, modelMl, modelPl, modelCl
            )
        err = sum(corrErr(world[:gradX], world[:tX]) +
        corrErr(world[:gradY], world[:tY]))
        println(i, " --- ",  err)
        world[:err] = err
        if err < 1E-6
            world[:itr] = i
            break
        end
    end
    return world
end

function produceSim(world)
    name = savename(world)
    println(name)
    res = runSim(world)
    # safesave(joinpath("..", "data", savename(world, "bson")), res)
    # println(res[:err], " --- ", name)
    return res
end

@time begin
    println("TASK: ", ENV["SLURM_ARRAY_TASK_ID"])

    # index = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
    # worldSet = dict_list(world)
    # worldSet = dict_list.(dict_list(world))
    # worldSet = collect(Iterators.flatten(worldSet))
    cosm = world
    cosm[:gain] = cosm[:ratio]/cosm[:stab]
    cosm[:loss] = (1-cosm[:ratio])/cosm[:stab]
    cosm[:fix] = string(cosm[:fixed][1], cosm[:fixed][2])

    w1 = produceSim(cosm)

    cosm[:nGens] = cosm[:realGen]
    # worldSet = dict_list(world)

    # for cosm in worldSet
    resWorld = produceSim(cosm)
    save(joinpath("/home", "mmp38", "rds", "hpc-work", savename(cosm, "bson")), resWorld)
end

# world = produceSim(worldSet[1])

# g = [0.2 -0.3 -0.11; -0.5 0.7 0.13]
# t = [0 0 0.5; 1 1 0.5]
# s = sign.(g)

# g[((t.==0) .& (s .> 0)) .| ((t.==1) .& (s .< 0)) .| (t .∉ [0,1])]



# for cosm in worldSet
#     produceSim(cosm)
# end

# for i in 1:world[:nGens]
#     println(i)
#     global world = updateMPC(copy(world), modelM, modelP, modelC, modelMl, modelPl, modelCl)

#     fFun(Fx, x) = callFun(
#         reshape(Fx, (world[:q], world[:n])), 
#         reshape(x, (world[:q], world[:n])), 
#         world[:tTr], 
#         world[:M], 
#         world[:P], 
#         world[:C], 
#         world[:d], 
#         world[:epsilon]
#         )

#     rFun(Fx, x) = callFunR(
#         reshape(Fx, (world[:q], world[:n])), 
#         reshape(x, (world[:q], world[:n])), 
#         solF.zero,
#         world[:tTr], 
#         world[:M], 
#         world[:P], 
#         world[:C], 
#         world[:d], 
#         world[:epsilon]
#         )

#     wFun(Fx, x) = callFunW(
#         reshape(Fx, (world[:q], world[:n])), 
#         reshape(x, (world[:q], world[:n])), 
#         solF.zero, 
#         world[:tTr], 
#         world[:Mf], 
#         world[:Ml], 
#         world[:Pf], 
#         world[:Pl],
#         world[:P],
#         world[:C],
#         world[:Cf],
#         world[:Cl],
#         world[:d], 
#         world[:epsilon]
#         )

    
#     solF = genSolF(fFun, world)
#     world[:tF] = solF.zero

#     solW = genSolW(wFun, world)
#     world[:tW] = solW.zero

#     solR = genSolR(rFun, world)
#     world[:tR] = solR.zero

#     gradX = substitute.(grads[1], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

#     gradY = substitute.(grads[2], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

#     world[:gradX] = Symbolics.value.(gradX)
#     world[:gradY] = Symbolics.value.(gradY)

#     world[:tX] = clamp.(
#         world[:tX] .+ world[:force]*world[:gradX],
#         # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
#         0,
#         1
#     )

#     world[:tY] = clamp.(
#         world[:tY] .+ world[:force]*world[:gradY],
#         # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
#         0,
#         1
#     )

#     display(world[:tF])
# end

# gradsMf = Symbolics.derivative.(wSysSelec, Mf)
# gradsPf = Symbolics.derivative.(wSysSelec, Pf)
# gradsCf = Symbolics.derivative.(wSysSelec, Cf)
# gradsMl = Symbolics.derivative.(wSysSelec, Ml)
# gradsPl = Symbolics.derivative.(wSysSelec, Pl)
# gradsCl = Symbolics.derivative.(wSysSelec, Cl)

# gradsF = [gradsMf, gradsPf, gradsCf]
# gradsFEq = [substitute.(x, atEqMPC) for x in gradsF]

# gradsL = [gradsMl, gradsPl, gradsCl]
# gradsLEq = [substitute.(x, atEqMPC) for x in gradsL]


# tF = fill!(zeros(world[:q], world[:n]), 1/world[:size])
# tM = fill!(zeros(world[:q], world[:n]), 0.1)
# tP = vcat([repeat([world[:k]+world[:b]*(i-1)], world[:n]) for i in 1:world[:q]]'...)
# world[:tTr] = fill!(zeros(world[:q], world[:q]), 0.1)
# tTr = fill!(zeros(world[:q], world[:q]), 0.1)
# tC = fill!(zeros(world[:q], world[:n]), 0.5)
# tX = fill!(zeros(world[:q], world[:n]), 0.5)
# tY = fill!(zeros(world[:q], world[:n]), 0.5)

# using NLsolve

# fFun(Fx, x) = callFun(
#     reshape(Fx, (world[:q], world[:n])), 
#     reshape(x, (world[:q], world[:n])), 
#     tTr, 
#     tM, 
#     tP, 
#     tC, 
#     world[:d], 
#     world[:epsilon]
#     )

# rFun(Fx, x) = callFunR(
#     reshape(Fx, (world[:q], world[:n])), 
#     reshape(x, (world[:q], world[:n])), 
#     solF.zero,
#     tTr, 
#     tM, 
#     tP, 
#     tC, 
#     world[:d], 
#     world[:epsilon]
#     )

# wFun(Fx, x) = callFunW(
#     reshape(Fx, (world[:q], world[:n])), 
#     reshape(x, (world[:q], world[:n])), 
#     solF.zero, 
#     tTr, 
#     tM, 
#     tM, 
#     tP,
#     tP,
#     tP, 
#     tC, 
#     tC, 
#     tC, 
#     world[:d], 
#     world[:epsilon]
#     )

# function makeFun(world) 
#     function f!(Fx, x)
#         ans = reshape(fFun(x), (world[:size], 1))
#         for i in 1:world[:size]
#             Fx[i] = ans[i]
#         end
#     end
#     return f!
# end
# solFun = makeFun(world)