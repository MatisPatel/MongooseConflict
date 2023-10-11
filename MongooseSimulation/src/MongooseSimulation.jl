module MongooseSimulation

export runSim, produceSim, produceOnceSim

# using ModelingToolkit 
using NLsolve
using Symbolics
using DrWatson
using BSON
using Statistics

#for testing ONLY 
# ENV["SLURM_ARRAY_TASK_ID"] = 10
# world = load(joinpath("tempDicts", string(ENV["SLURM_ARRAY_TASK_ID"], ".bson")))

# world = Dict{Symbol, Any}(
#     :nGens => 1,
#     :realGen => 3,
#     :worldSize => [3, 3],
#     :ratio => 0.5,
#     :stab => 5,
#     :fixed => [1,2],
#     :gain => 0.1,
#     :loss => 0.1,
#     :basem => 0.1,
#     :k => 0.5,
#     :b => 1.0,
#     :d => 0.5,
#     :epsilon => 5,
#     :multX => 0.1,
#     :multY => 0.1,
#     :shape_X_cost => 0.5, 
#     :shape_Y_cost => 0.5,
#     :grad_rate => 0.1,
#     :learning_rate => 0.01,
#     :decay => 0.9,
#     :SolFFails => 0,
#     :SolWFails => 0,
#     :SolRFails => 0,
#     :saveKeys => [
#         :multX,
#         :multY,
#         :shape_X_cost,
#         :shape_Y_cost,
#         :b,
#         :k,
#         :d,
#         :q,
#         :n,
#         :decay, 
#         :epsilon, 
#         :learning_rate, 
#         :ratio, 
#         :stab,]
# )

# @variables g, l, d, epsilon, k, b, B, deltaX, deltaY, dummy, t

# # world = load(joinpath("/home", "mmp38", "mongooseConflict", "code", "tempDicts", string(ENV["SLURM_ARRAY_TASK_ID"], ".bson")))
# world[:q] = world[:worldSize][1] 
# world[:n] = world[:worldSize][2] 
# world[:size] = world[:n]*world[:q] 
# world[:gradX] = zeros(world[:q], world[:n])
# world[:gradY] = zeros(world[:q], world[:n])
# world[:step_size_x] = zeros(world[:q], world[:n])
# world[:step_size_y] = zeros(world[:q], world[:n])
# world[:err_list] = zeros(world[:nGens])

# @variables begin 
#     F[1:world[:q], 1:world[:n]]
#     W[1:world[:q], 1:world[:n]]
#     R[1:world[:q], 1:world[:n]]
#     M[1:world[:q], 1:world[:n]]
#     Mf[1:world[:q], 1:world[:n]]
#     Ml[1:world[:q], 1:world[:n]]
#     P[1:world[:q], 1:world[:n]]
#     Pf[1:world[:q], 1:world[:n]]
#     Pl[1:world[:q], 1:world[:n]]
#     C[1:world[:q], 1:world[:n]]
#     Cf[1:world[:q], 1:world[:n]]
#     Cl[1:world[:q], 1:world[:n]]
#     Tr[1:world[:q], 1:world[:q]]
#     X[1:world[:q], 1:world[:n]]
#     Y[1:world[:q], 1:world[:n]]
#     Xf[1:world[:q], 1:world[:n]]
#     Yf[1:world[:q], 1:world[:n]]
#     Xl[1:world[:q], 1:world[:n]]
#     Yl[1:world[:q], 1:world[:n]]
# end

function try_function(fun)
    attempts = 0
    while attempts < 10
        attempts += 1
        try
            result = fun()
            return result
        catch
            println("Attempt $attempts failed.")
        end
    end
    println(string("Function ", fun, " failed after 10 attempts."))
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
                if newQ == q-1
                    Fp[q, n] -= world[:loss] * F[q, n]
                    Fp[newQ, n] += world[:loss] * F[q, n]
                elseif newQ == q+1 
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
                if newQ == q-1
                    Wn[q, n] += world[:loss] * W[newQ, n]
                    Wd[q, n] += world[:loss]
                elseif newQ == q+1  
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
                if newQ == q-1
                    Rp[q, n] -= world[:loss] * F[q, n] * R[q, n]
                    Rp[newQ, n] += world[:loss] * F[q, n] * R[q, n]
                elseif newQ == q+1 
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
    pbar = sum( F.*P * [n-1 for n in 1:world[:n]])
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
            Rp[q, n] -= d * Pbar * F[q, n] * R[q, n]
            Rp[q, n+1] += d * Pbar * F[q, n] * ((newN-2)/newN)*R[q, n]            
        end
    end
    return Rp
end

function victory(a, b)
    return (a+1E-8)*(a+b+2*1E-8)^-1.0
end

function addFights(Fp, F, C, epsilon, world)
    for q in 1:world[:q]
        for n in 1:world[:n]
            for qOpp in 1:world[:q]
                for nOpp in 1:world[:n]
                    if !(n==1 && nOpp==1)
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
                    else
                        # println("NO MATCH") 
                    end
                end
                end
            end
        end
    end                                    
    return Fp
end

function addFights(Wn, Wd, W, F, C, Cf, Cl, epsilon, world)
    for q in 1:world[:q]-1
        for n in 2:world[:n]
            for qOpp in 2:world[:q]
                for nOpp in 1:world[:n]
                    # println("fight", F[q,n], " -> ", F[qOpp, nOpp])
                    pEncounter = F[qOpp, nOpp] * epsilon
                    pVictory = victory(
                        Cf[q, n] + (n-2) * Cl[q, n], 
                        (nOpp-1) * C[qOpp, nOpp]
                    )
                    # pVictory = victory(0.1*(n-1), 0.1*(nOpp-1))
                    # println(pVictory)
                    fWin = pEncounter * pVictory
                    fLoss = pEncounter * (1 - pVictory)
                    Wn[q, n] +=  fWin * W[q + 1, n]
                    Wd[q, n] += fWin 
                end 
            end 
        end 
    end
    for q in 2:world[:q]
        for n in 1:world[:n]
            for qOpp in 1:world[:q]-1
                for nOpp in 2:world[:n]
                    # println("fight", F[q,n], " -> ", F[qOpp, nOpp])
                    pEncounter = F[qOpp, nOpp] * epsilon
                    pVictory = victory( 
                    (nOpp-1) * C[qOpp, nOpp],
                    Cf[q, n] + (n-2) * Cl[q, n]
                    )
                    # println(pVictory)
                    # pVictory = victory(0.1*(nOpp-1), 0.1*(n-1))
                    fWin = pEncounter * pVictory
                    fLoss = pEncounter * (1 - pVictory)
                    Wn[q, n] +=  fWin * W[q - 1, n]
                    Wd[q, n] += fWin 
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
# function addFightsR(Rp, R, F, C, epsilon, world)
#     for q in 1:world[:q]
#         for n in 1:world[:n]
#             for qOpp in 1:world[:q]
#                 for nOpp in 1:world[:n]
#                     # println("fight", F[q,n], " -> ", F[qOpp, nOpp])
#                     # halved to not double count x->y and y->x
#                     pEncounter  =  F[q, n] * F[qOpp, nOpp] * epsilon
#                     pVictory    =  victory((n-1) * C[q, n], (nOpp-1) * C[qOpp, nOpp])
#                     # print(pVictory)
#                     fWin = pEncounter*pVictory
#                     fLoss = pEncounter*(1-pVictory)
#                     # fight logic 
#                     # what is focal patch is minimal richness? then can only win losses dont change freq
#                     # or if qOpp is maximum and q is less than it
#                     if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
#                         # println("only win")
#                         Rp[q,n] -= fWin * R[q, n]
#                         Rp[q+1,n] += fWin * R[q, n]
#                     # what if focal patch is non-zero but other patch is 0. Then can only lose
#                     elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
#                         # println("only lose")
#                         Rp[q,n] -= fLoss * R[q,n]
#                         Rp[q-1,n] += fLoss * R[q,n]
#                     # if equal richness and neither max not min then they can win or lose 
#                     elseif (q == qOpp) & (q > 1) & (q < world[:q])
#                         # println("both win or lose")
#                         # q wins
#                         Rp[q,n] -= fWin * R[q,n]
#                         Rp[q+1,n] += fWin * R[q,n]
#                         # q loses
#                         Rp[q,n] -= fLoss * R[q,n]
#                         Rp[q-1,n] += fLoss * R[q,n]
#                     # or if diffrent and not min or max
#                     elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
#                         # println("both win or lose")
#                         # q wins
#                         Rp[q,n] -= fWin * R[q,n]
#                         Rp[q+1,n] += fWin * R[q,n]
#                         # q loses
#                         Rp[q,n] -= fLoss * R[q,n]
#                         Rp[q-1,n] += fLoss * R[q,n]
#                     # else
#                         # println("NO MATCH") 
#                     end
#                 end
#             end
#         end
#     end                                    
#     return Rp
# end

function makeFsys(F, M, P, C, d, world)
    Fp = makeFArray(world)
    Fp1 = addTrans(copy(Fp), F, world)
    Fp2 = addMortality(copy(Fp1), F, M, world)
    Fp3 = addLocBirth(copy(Fp2), F, P, d, world)
    Fp4 = addImmigration(copy(Fp3), F, P, d, world)
    Fp5 = addFights(copy(Fp4), F, C, epsilon, world)
    Fp5[1,1] = 1 - sum(F) + F[1,1]
    return Fp5
end

function makeFocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    println(size(model))
    println(world[:q], " ", world[:n])
    for q in 1:world[:q]
        for n in 2:world[:n]
            println("q: ", q, " n: ", n)
            # model[q, n] = world[:basem] * 
            # exp(-1 * ((Xl[q,n]*(n-2) + 
            # Xf[q,n]))) + 
            # world[:multX]*Xf[q,n]^world[:shape_X_cost] + 
            # world[:multY]*Yf[q,n]^world[:shape_Y_cost]
        end
    end
    return model
end

function makeLocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = world[:basem] * exp(-1 * ((Xl[q,n]*(n-2) + Xf[q,n]))) + world[:multX]*Xl[q,n]^world[:shape_X_cost] + world[:multY]*Yl[q,n]^world[:shape_Y_cost]
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
    end
    return model
end

function makeLocalModelP(world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = prodFun(q, n, world)
        end
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
        for n in 2:world[:n]
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
    wSys = Symbolics.scalarize(Wn6./Wd6.-W)
    wSys[1,2] = 1-W[1,2]
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

function makeModelExpr(world)
    modelMf = makeFocalModelM(Xf, Xl, Yf, Yl, world)
    modelMl = makeLocalModelM(Xf, Xl, Yf, Yl, world)

    modelPf = makeFocalModelP(world)
    modelPl = makeLocalModelP(world)

    modelCf = makeFocalModelC(Yf, world)
    modelCl = makeLocalModelC(Yl, world)

    atEqXY = Dict([Xf => X, Xl => X, Yf => Y, Yl => Y])
    modelM = substitute(modelMf, atEqXY)
    modelP = substitute(modelPf, atEqXY)
    modelC = substitute(modelCf, atEqXY)

    # modelMFun = eval(ModelingToolkit.build_function(
    #     modelM, B, X, Y
    # )[1])
    # modelPFun = eval(ModelingToolkit.build_function(
    #     modelP, b, k
    # )[1])
    # modelCFun = eval(ModelingToolkit.build_function(
    #     modelC, Y
    # )[1])

    # modelMlFun = eval(ModelingToolkit.build_function(
    #     modelM, B, Xf, Xl, Yl
    # )[1])
    # modelPlFun = eval(ModelingToolkit.build_function(
    #     modelP, b, k
    # )[1])
    # modelClFun = eval(ModelingToolkit.build_function(
    #     modelC, Yl
    # )[1])
    return modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC
end

function updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)
    world[:M] = Symbolics.value.(substitute(modelM, Dict([B=>world[:basem], X => world[:tX], Y => world[:tY]]))) 
    world[:P] = Symbolics.value.(modelP)
    world[:C] = Symbolics.value.(substitute(modelC, Dict([Y => world[:tY]])))
    world[:Mf] = world[:M]
    world[:Pf] = world[:P]
    world[:Cf] = world[:C]
    world[:Ml] = Symbolics.value.(substitute(modelMl, 
        Dict([
            B => world[:basem], 
            Xf => world[:tX], 
            Xl => world[:tX], 
            Yl => world[:tY],  
        ])
    ))
    world[:Pl] = Symbolics.value.(modelPl)
    world[:Cl] = Symbolics.value.(substitute(modelCl, Dict([Yl => world[:tY]])))
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
    wSysXY = substitute(wSysSelec, Dict([
        Ml =>modelMl,
        Mf =>modelMf, 
        P =>modelP, 
        Pf =>modelPf, 
        Pl =>modelPl,
        C =>modelC, 
        Cf =>modelCf, 
        Cl =>modelCl, 
    ])
    );
    wDiffXf = Symbolics.derivative.(wSysXY, Xf);
    wDiffXl = Symbolics.derivative.(wSysXY, Xl);
    wDiffYf = Symbolics.derivative.(wSysXY, Yf);
    wDiffYl = Symbolics.derivative.(wSysXY, Yl);
    inclusiveX = wDiffXf .+ R .* wDiffXl
    inclusiveY = wDiffYf .+ R .* wDiffYl

    Eq = Dict([
        Xf =>X, 
        Xl =>X, 
        Yf =>Y, 
        Yl =>Y,
        d => world[:d],
        epsilon  => world[:epsilon]
    ]
    )

    incXEq = substitute(Symbolics.scalarize(inclusiveX), Eq);

    incYEq = substitute(Symbolics.scalarize(inclusiveY), Eq);

    grads = [incXEq, incYEq]
    directSel = [wDiffXf, wDiffYf]
    indirectSel = [R .* wDiffXl, R .* wDiffYl]

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
        transpose(
            reshape(
                repeat(
                    vcat([0, 0], [1/p for p in 2:(world[:n]-1)]), 
                    world[:q]
                ), 
                (world[:q], world[:n])
            )
        )
    )
end

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
        world[:Cl],
        world[:Cf],
        world[:d], 
        world[:epsilon]
        )

    solF = zeros(world[:q], world[:n])
    try
        try 
            solF = genSolF2(fFun, world)
        catch
            solF = genSolF(fFun, world)
        finally 
            # println("F sol failed for:\n", 
            #         savename(world))
            # display(world[:tF])
            # display(world[:tW])
            
        end
    catch
        # println("Defaulting to previous F")
        # display(world[:tF])
        solF = world[:tF]
    end
    # solF = genSolF(fFun, world)
    # println(solF.iterations)
    # display(solF.zero)
    if hasproperty(solF, :zero)
        world[:tF] = solF.zero;
    else 
        world[:tF] = world[:tF]
        world[:SolFFails] += 1
    end
    
    solW = zeros(world[:q], world[:n])
    try
        try 
            solW = genSolW2(wFun, world)
        catch
            solW = genSolW(wFun, world)
        finally 
            # println("W sol failed for:\n", 
            #         savename(world))
            # display(world[:tF])
            # display(world[:tW])
            
        end
    catch
        # println("Defaulting to previous W")
        # display(world[:tW])
        solW = world[:tW]
    end
    # solW = genSolW(wFun, world)
    # println(solW.iterations)
    # display(solW.zero)
    if hasproperty(solW, :zero)
        world[:tW] = solW.zero;
    else 
        world[:tW] = world[:tW]
        world[:SolWFails] += 1
    end

    solR = zeros(world[:q], world[:n])
    try
        try 
            solR = genSolR2(rFun, world)
        catch
            solR = genSolR(rFun, world)
        finally
            # println("R sol failed for:\n", savename(world))
            # display(world[:tF])
            # display(world[:tW])
        end
    catch 
        # println("Defaulting to previous R")
        # display(world[:tR])
        solR = world[:tR]
    end
    
    # solR = genSolR(rFun, world)
    # println(solR.iterations)
    # display(solR.zero)
    if hasproperty(solR, :zero)
        world[:tR] = solR.zero;
    else 
        world[:tR] = world[:tR]
        world[:SolRFails] += 1
    end

    # attempt to solve for F 10 times 
    # attempts = 0 
    # solF = :missing
    # while attempts < 10
    #     try
    #         solF = genSolF(fFun, world)
    #         break
    #     catch
    #         attempts += 1
    #     end
    # end
    # if solF != :missing
    #     world[:tF] = solF.zero
    # else
    #     return world, false
    # end

    # attempt to solve for W 10 times 
    # attempts = 0    
    # solW = :missing
    # while attempts < 10
    #     try
    #         solW = genSolW(wFun, world)
    #         break
    #     catch
    #         attempts += 1
    #     end
    # end
    # if solW != :missing
    #     world[:tW] = solW.zero
    # else
    #     return world, false
    # end
    # solW = genSolW(wFun, world)
    # world[:tW] = solW.zero

    # attempt to solve for R 10 times
    # attempts = 0
    # solR = :missing
    # while attempts < 10
    #     try
    #         solR = genSolR(rFun, world)
    #         break
    #     catch
    #         attempts += 1
    #     end
    # end
    # if solR != :missing
    #     world[:tR] = clamp.(solR.zero, 0, 1)
    # else
    #     return world, false
    # end 

    # solR = zeros(world[:q], world[:n])
    # try 
    #     solR = genSolR2(rFun, world)
    # catch
    # solR = genSolR(rFun, world)
    # # end 
    # world[:tR] = clamp.(solR.zero, 0, 1)

    gradX = substitute(grads[1], Dict([X =>world[:tX], Y =>world[:tY], F =>world[:tF], W =>world[:tW], R =>world[:tR]]))

    gradY = substitute(grads[2], Dict([X =>world[:tX], Y =>world[:tY], F =>world[:tF], W =>world[:tW], R =>world[:tR]]))

    # world[:tgradX] = world[:gradX]
    # world[:tgradY] = world[:gradY]
    world[:gradX] = clamp.(Symbolics.value.(gradX), -1, 1)
    world[:gradY] = clamp.(Symbolics.value.(gradY), -1, 1)

    # (world[:step_size_x], world[:signs_x]) = update_step_size(world[:tgradX], world[:gradX], world[:step_size_x], world)
    # (world[:step_size_y], world[:signs_y]) = update_step_size(world[:tgradY], world[:gradY], world[:step_size_y], world)

    # world[:tX] = clamp.(world[:tX], 0, 1) 
    # world[:tX] = clamp.(world[:tX] .+ world[:step_size_max] .* world[:signs_x], 0, 1)
    # world[:tX] = clamp.(
    #     world[:tX] .+ world[:step_size_max]*clamp.(world[:gradX], -1, 1),
    #     # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
    #     0,
    #     1
    # )


    if world[:err] >= 1 && !world[:randStart]

    world[:tX] = clamp.(
        world[:tX] .+ world[:grad_rate]*clamp.(world[:gradX], -1, 1),
        # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
        1E-8,
        1
    )

    world[:tY] = clamp.(
        world[:tY] .+ world[:grad_rate]*clamp.(world[:gradY], -1, 1),
        # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
        1E-8,
        1
    )

    else
        (world[:update_x], world[:cache_x]) = RMSprop_update(world[:learning_rate], world[:cache_x], world[:cache_y], world[:gradX], world[:decay])

        (world[:update_y], world[:cache_y]) = RMSprop_update(world[:learning_rate], world[:cache_y], world[:cache_x], world[:gradY], world[:decay])

        world[:tX] = clamp.(world[:tX] .+ world[:update_x], 1E-8, 1)
        world[:tY] = clamp.(world[:tY] .+ world[:update_y], 1E-8, 1)
    end

    return world, false
end

function runSim(world)
    world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
    world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))

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

    for i in 1:world[:nGens]
        world[:itr] = i
        (world, e_flag) = step(world, callFun, callFunW, callFunR,
            grads,
            modelM, modelP, modelC, modelMl, modelPl, modelCl
            )
        if e_flag == true 
            println(savename(world), " \nSTEP HAD A FATAL ERROR");
            break
        end

        err = sum(corrErr(world[:gradX], world[:tX]) +
        corrErr(world[:gradY], world[:tY]))
        if world[:verbose]
            println(i, " --- ",  err)
        end
        world[:err] = err
        world[:err_list][i] = err
        if err < 1E-6
            world[:itr] = i
            break
        end
        # if i%1000 == 0
        #     save(joinpath("/home", "mmp38", "rds", "hpc-work", savename(cosm, "bson")), world)
        # end
    end
    return world
end

function testRatios!(resWorld)
    for newRatio in 0.1:0.05:0.9
        newWorld = deepcopy(resWorld)
        newWorld[:ratio] = newRatio 
        newWorld[:force] = 0.0 
        newWorld[:gain] = newWorld[:ratio]/newWorld[:stab]
        newWorld[:loss] = (1-newWorld[:ratio])/newWorld[:stab]
        newWorld[:nGens] = 1 
        newWorld[:tF] = reshape(
            repeat([1/newWorld[:size]], newWorld[:size]), (newWorld[:q], newWorld[:n])
        )
        newWorld[:tR] = reshape(
            repeat(
                hcat([0 0], [1/p for p in 2:(newWorld[:n]-1)]'), 
                newWorld[:q]
            ), 
            (newWorld[:q], newWorld[:n])
        )

        modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC = makeModelExpr(newWorld)

        newWorld = updateMPC(newWorld, modelM, modelP, modelC, modelMl, modelPl, modelCl)

        # create F system
        fSys = makeFsys(F, M, P, C, d, newWorld)
        funF = build_function(
            fSys, F, M, P, C, d, epsilon;
            expression=Val{false}
            );
        callFun = eval(funF[2]);
        newWorld[:itr] = -1
        fFun(Fx, x) = callFun(
            reshape(Fx, (newWorld[:q], newWorld[:n])), 
            reshape(x, (newWorld[:q], newWorld[:n])),  
            newWorld[:M], 
            newWorld[:P], 
            newWorld[:C], 
            newWorld[:d], 
            newWorld[:epsilon]
            )
        solF = genSolF(fFun, newWorld)
        resWorld[Symbol(replace(string("rF","_", newWorld[:ratio]), "."=>"_"))] = solF.zero
    end
    return resWorld
end

function produceSim(world, save=false)
    name = savename(world)
    if world[:verbose]
        println(name)
    end
    res = runSim(world)
    if save == true
        wsave(joinpath("..", "data", savename(world, "bson")), res)
    end
    println(res[:err], " --- ", name)
    return res
end

function produceOnceSim(world, save=false)

    name = savename(world, "bson", accesses=world[:saveKeys])
    if name in readdir(joinpath("..", "data"))
        println("skipping ", name)
        return world
    end
    (res, elapsed_time) = @timed runSim(world) 
    # (res, error_flag) = result
    println(string("Time: ", elapsed_time))
    res[:timing] = elapsed_time
    if save == true 
        if world[:verbose]
            println(string("Saving sim ", name))
            println(res[:itr], " --- ", res[:err])
        end
        safesave(joinpath("..", "data", name), res);
    end
    
    return res
end

# @time begin
#     # println("TASK: ", ENV["SLURM_ARRAY_TASK_ID"])

#     # index = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
#     # worldSet = dict_list(world)
#     # worldSet = dict_list.(dict_list(world))
#     # worldSet = collect(Iterators.flatten(worldSet))
#     cosm = world
#     cosm[:gain] = cosm[:ratio]/cosm[:stab]
#     cosm[:loss] = (1-cosm[:ratio])/cosm[:stab]
#     cosm[:fix] = string(cosm[:fixed][1], cosm[:fixed][2])

#     w1 = produceSim(cosm)s

#     cosm[:nGens] = cosm[:realGen]
#     # worldSet = dict_list(world)

#     # for cosm in worldSet
#     resWorld = produceSim(cosm)
#     finalWorld = testRatios!(resWorld)
#     # save(joinpath("/home", "mmp38", "rds", "hpc-work", savename(cosm, "bson")), finalWorld)
# end

function update_step_size(prev_grad, curr_grad, step_size, world)
    # println(prev_grad, curr_grad)
    signs = sign.(curr_grad .* prev_grad)
    for i in 2:world[:q]
        for j in 2:world[:n] 
            if signs[i, j] > 0 
                step_size[i, j] = min(step_size[i,j] * world[:inc_factor], world[:step_size_max]) 
            elseif signs[i, j] < 0 
                step_size[i, j] = max(step_size[i,j] * world[:dec_factor], world:[:step_size_min])
            end
        end
    end
    return step_size, signs
end

function RMSprop_update(learning_rate, cache, grad, decay)
    cache = decay .* cache .+ (1 - decay) .* grad.^2
    return grad .* (learning_rate ./ ((sqrt.(cache)) .+ 1E-8)), cache
end

function RMSprop_update(learning_rate, cache1, cache2, grad, decay)
    cache1 = decay .* cache1 .+ (1 - decay) .* clamp.(grad, -1, 1).^2
    return grad .* (learning_rate ./ ((sqrt.(0.5.*(cache1 .+ mean(cache2)))) .+ 1E-8)), cache1
end

end