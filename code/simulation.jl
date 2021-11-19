using ModelingToolkit 
using NLsolve
using Symbolics
using DrWatson

@variables g, l, d, epsilon, k, b, B, deltaX, deltaY, dummy, t

pprint(x) = round.(x, digits=4)

world = Dict{Symbol, Any}(
    :force => [0.0],
    :itr =>0,
    :nGens => 1,
    :realGen =>  [@onlyif(:force != 0, 10000), @onlyif(:force==0, 1)],
    :q => 2,
    :n => 2,
    # :gain => [0.05, 0.1, 0.15, 0.2, 0.25],
    # :loss => [0.05, 0.1, 0.15, 0.2, 0.25],
    :stab => [5],
    # :ratio => collect(0.1:0.1:0.9),
    :ratio => [0.5],
    :basem => 0.1,
    :k => 0.1,
    :b => 0.3,
    :d => [0.5],
    :epsilon => 5,
    :multX => 0.1,
    :multY => 0.1,
    :fixed => [[1,2]],
    # :fixed => [[1,2], [5,3], [3,3], [3,2]],
    :placeholder => [[], []],
    :solver => :normal, #normal, ADAM, momentum, weighted
    :momentum => 0.9,
    :beta1 => 0.8,
    :beta2 => 0.999
)

world[:size] = world[:n]*world[:q]

save("tempDicts/test.bson", world)

worldSet = dict_list(world)

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

# function addDistantBirth(Wn, Wd, W, F, Pf, world)
#     for q in 1:world[:q]
#         for n in 2:world[:n]
#             for qNew in 1:world[:q]
#                 for nNew in 1:(world[:n]-1)
#                     Wn[q, n] += d * Pf[q,n] * F[qNew, nNew] * (W[q,n] + W[qNew, nNew+1])
#                     Wd[q, n] += d * Pf[q,n] * F[qNew, nNew]
#                 end
#             end
#         end
#     end
#     return Wn, Wd
# end

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

# function addDistantBirth(Wn, Wd, W, F, Pf, world)
#     for q in 1:world[:q]
#         for n in 2:world[:n]
#             for qNew in 1:world[:q]
#                 for nNew in 1:(world[:n])
#                     if nNew != world[:n]
#                         Wn[q, n] += d * Pf[q,n] * F[qNew, nNew] * (W[q,n] + W[qNew, nNew+1])
#                         Wd[q, n] += d * Pf[q,n] * F[qNew, nNew]
#                     else 
#                         W[q, n] += d * Pf[q,n] * F[qNew, nNew]*W[q,n]
#                         Wd[q,n] += d * Pf[q,n] * F[qNew, nNew]
#                     end
#                 end
#             end
#         end
#     end
#     return Wn, Wd
# end

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
                    
                    # if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
                    #     # println("only win")
                    #     Wn[q, n] += fWin * W[q + 1, n]
                    #     Wd[q, n] += fWin
                    # # what if focal patch is non-zero but other patch is 0. Then can only lose
                    # elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
                    #     # println("only lose")
                    #     Wn[q, n] += fLoss * W[q - 1, n]
                    #     Wd[q, n] += fLoss
                    # # if equal richness and neither max not min then they can win or lose 
                    # elseif (q == qOpp) & (q > 1) & (q < world[:q])
                    #     # println("both win or lose")
                    #     # q wins
                    #     Wn[q, n] += fWin * W[q + 1, n]
                    #     Wd[q, n] += fWin 
                    #     # q loses
                    #     Wn[q, n] += fLoss * W[q - 1, n]
                    #     Wd[q, n] += fLoss
                    # # or if diffrent and not min or max
                    # elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
                    #     # println("both win or lose")
                    #     Wn[q, n] += fWin * W[q + 1, n]
                    #     Wd[q, n] += fWin 
                    #     # q loses
                    #     Wn[q, n] += fLoss * W[q - 1, n]
                    #     Wd[q, n] += fLoss
                    # else
                    #     println("NO MATCH") 
                    # end
    #             end
    #         end
    #     end
    # end
#     return Wn, Wd
# end

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
    # Fp5 =Fp4
    Fp5 = addFights(copy(Fp4), F, C, epsilon, world)
    Fp5[1,2] = 1 - sum(F)
    return Fp5
end

function mortFun2(n, xf, xl, yf, yl, world)
    calc = B * 2.718^(-1 * (n-1)*xl) + world[:multX]*xf^2 + world[:multY]*yf^2
    return calc
end 

# function mortFun(n, xf, xl, yf, yl, world)
#     calc = world[:basem] * exp(-1 * ((xl*(n-1) + xf))) + world[:multX]*xf^2 + world[:multY]*yf^2
#     return calc
# end 

function mortFun(n, xf, xl, yf, yl, world)
    calc = world[:basem] * exp(-1 * ((xl*(n-1) + xf))) + world[:multX]*xf^2 + world[:multY]*yf^2
    return calc
end 

function makeFocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            # model[q, n] = mortFun(n-1, Xf[q, n], Xl[q, n], Yf[q, n], Yl[q, n], world)
            model[q, n] = world[:basem] * exp(-1 * ((Xl[q,n]*(n-2) + Xf[q,n]))) + world[:multX]*Xf[q,n]^2 + world[:multY]*Yf[q,n]^2
        end
    end
    return model
end
# makeFocalModelM(world[:tX], world[:tX], world[:tX], world[:tX], world)

function makeLocalModelM(Xf, Xl, Yf, Yl, world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            # model[q, n] = mortFun(n-1, Xl[q, n], Xf[q, n], Yl[q, n], Yf[q, n], world)
            model[q, n] = world[:basem] * exp(-1 * ((Xl[q,n]*(n-2) + Xf[q,n]))) + world[:multX]*Xl[q,n]^2 + world[:multY]*Yl[q,n]^2
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
        # model[q, 1] = 0.0
    end
    return model
end

function makeLocalModelP(world)
    model = makeFArray(world)
    for q in 1:world[:q]
        for n in 2:world[:n]
            model[q, n] = prodFun(q, n, world)
        end
        # model[q, 1] = 0.0
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
    wSys = (Wn6./Wd6).-W
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
    # atEqXY = Dict([Xf .=> X, Xl .=> X, Yf .=> Y, Yl .=> Y])
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

function makeSelGrads2(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, 
    modelPf, modelPl, modelCf, modelCl, world)
    wSysXY = substitute.(wSysSelec, matSub(
        Mf.=>modelMf, 
        Ml.=>modelMl,
        P.=>modelP, 
        Pf.=>modelPf, 
        Pl.=>modelPl,
        C.=>modelC, 
        Cf.=>modelCf, 
        Cl.=>modelCl, 
        ))
    wDiffXf = Symbolics.derivative.(wSysXY, Xf)
    wDiffXl = Symbolics.derivative.(wSysXY, Xl)
    wDiffYf = Symbolics.derivative.(wSysXY, Yf)
    wDiffYl = Symbolics.derivative.(wSysXY, Yl)
    inclusiveX = wDiffXf .+ R .* wDiffXl
    inclusiveY = wDiffYf .+ R .* wDiffYl

    Eq = matSub(
        Xf.=>X, 
        Xl.=>X, 
        Yf.=>Y, 
        Yl.=>Y,
        [d].=>[world[:d]],
        [epsilon].=>[world[:epsilon]]
    )

    incXEq = substitute.(inclusiveX, Eq)

    incYEq = substitute.(inclusiveY, Eq)

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
            reshape(
                repeat(
                    vcat([0, 0], [1/p for p in 2:(world[:n]-1)]), 
                    world[:q]
                ), 
                (world[:q], world[:n])
            )
        )
end

# function genSolR(rFun, world)
#     return nlsolve(
#         rFun, ones(world[:q], world[:n]).*0.1
#         )
# end


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
    world[:itr] += 1
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

    
    solF = genSolF(fFun, world)
    world[:tF] = solF.zero

    solW = genSolW(wFun, world)
    world[:tW] = solW.zero
    world[:solW] = solW

    solR = genSolR(rFun, world)
    world[:tR] = solR.zero

    gradX = substitute.(grads[1], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

    gradY = substitute.(grads[2], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

    world[:gradX] = Symbolics.value.(gradX)
    world[:gradY] = Symbolics.value.(gradY)

    if world[:solver] == :momentum
        if !(haskey(world, :gradXLast))
            world[:gradXLast] = world[:gradX].*0.0
            world[:gradYLast] = world[:gradY].*0.0
            world[:gradXCurr] = world[:gradX]
            world[:gradYCurr] = world[:gradY]
        end 

        world[:gradXLast] = world[:gradYCurr]
        world[:gradYLast] = world[:gradYCurr]
        world[:gradXCurr] = world[:gradX] .+ world[:momentum] .* world[:gradXLast]
        world[:gradYCurr] = world[:gradY] .+ world[:momentum] .* world[:gradYLast]

        world[:tX] = clamp.(
            world[:tX] .+ world[:force]*clamp.(world[:gradXCurr], -1, 1),
            # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
            0,
            1
        )

        world[:tY] = clamp.(
            world[:tY] .+ world[:force]*clamp.(world[:gradYCurr], -1, 1),
            # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
            0,
            1
        )
    end

    if world[:solver] == :weighted
        world[:tX] = clamp.(
            world[:tX] .+ (0.01 .+ world[:tW]).*(world[:force]*world[:gradX]),
            # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
            0,
            1
        )

        world[:tY] = clamp.(
            world[:tY] .+ (0.01 .+ world[:tW]).*(world[:force]*world[:gradY]),
            # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
            0,
            1
        )
    end
    # 980 673
    if world[:solver] == :normal
        world[:tX] = clamp.(
            world[:tX] .+ (world[:force]*world[:gradX]),
            # clamp.(world[:force]*world[:gradX], -world[:force], world[:force]), 
            0,
            1
        )

        world[:tY] = clamp.(
            world[:tY] .+ (world[:force]*world[:gradY]),
            # clamp.(world[:force]*world[:gradY], -world[:force], world[:force]), 
            0,
            1
        )
    end

    if world[:solver] == :adam 
        world[:mX] = world[:beta1] .* world[:mX] .+ (1.0 - world[:beta1])*world[:gradX]
        world[:vX] = (world[:beta2] .* world[:vX]) .+ (1.0 - world[:beta2]).*world[:gradX].^2
        mXhat = world[:mX]./(1-world[:beta1]^world[:itr])
        vXhat = world[:vX]./(1-world[:beta2]^world[:itr])
        world[:tX] = clamp.(
            world[:tX] .+ (world[:force] .* mXhat ./ (sqrt.(vXhat).+1E-8)),
            0,
            1
        )

        world[:mY] = world[:beta1] .* world[:mY] .+ (1.0 - world[:beta1])*world[:gradY]
        world[:vY] = (world[:beta2] .* world[:vY]) .+ (1.0 - world[:beta2]).*world[:gradY].^2
        mYhat = world[:mY]./(1-world[:beta1])
        vYhat = world[:vY]./(1-world[:beta2])
        world[:tY] = clamp.(
            world[:tY] .+ (world[:force] .* mYhat ./ (sqrt.(vYhat).+1E-8)),
            0,
            1
        )
    end

    return world
end

function runSim(world)
    world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-10))
    world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-10))
    world[:tTr] = simpleTrans(world)

    if world[:solver] == :adam 
        world[:mX] = fill!(zeros(world[:q], world[:n]), 0.0)
        world[:vX] = fill!(zeros(world[:q], world[:n]), 0.0)
        world[:mY] = fill!(zeros(world[:q], world[:n]), 0.0)
        world[:vY] = fill!(zeros(world[:q], world[:n]), 0.0)
    end


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

    # Wavg = copy(W) 
    # for q in 1:world[:q]
    #     for n in 1:world[:n]
    #         Wavg[q,n] = 1 - (sum(W.*F) - W[q,n])*F[q,n]
    #     end
    # end

    # wSys = substitute.(wSys, matSub(W.=>Wavg))

    # wSys[world[:fixed]...] = 1-W[world[:fixed]...]
    for q in 1:world[:q]
        wSys[q, 1] = W[q,1]
    end

    wSys[world[:fixed]...] = 1 - sum((W.*F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))')./sum(F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))'))

    # wSys = substitute.(wSys, matSub(W[:,1].=>0.0))

    wSysSelec = Wn./Wd
    funW = ModelingToolkit.build_function(
        wSys, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
        expression=Val{false}
        );
    callFunW = eval(funW[2]);

    grads, directSel, indirectSel = makeSelGrads2(wSysSelec, modelM, modelP, modelC, modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, world);

    for i in 1:world[:nGens]
        world = step(world, callFun, callFunW, callFunR,
            grads,
            modelM, modelP, modelC, modelMl, modelPl, modelCl
            )
        err = sum(corrErr(world[:gradX], world[:tX]) +
        corrErr(world[:gradY], world[:tY]))
        println(i, " --- ",  err)
        # println("------ ", sum(world[:tF][:, 2:end]))
        # if err < 1E-6 
        #     world[:force] = world[:force]/10
        # end
        if err < 1E-8
            EqVals = matSub(X.=>world[:tX], Y.=>world[:tY], Xf.=>world[:tX], Yf.=>world[:tY], Xl.=>world[:tX], Yl.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR], epsilon.=>[world[:epsilon]], d.=>[world[:d]])
            world[:gradXdir] = substitute.(directSel[1], EqVals)
            world[:gradYdir] = substitute.(directSel[2], EqVals)
            world[:gradXind] = substitute.(indirectSel[1], EqVals)
            world[:gradYind] = substitute.(indirectSel[2], EqVals)
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


function calcFights(world)
    #number of valid fights 
    totalF = 0 
    partF = 0
    for q in 1:world[:q]
        for n in 1:world[:n]
            for qOpp in 1:world[:q]
                for nOpp in 1:world[:n]
                    if !(n==1 && nOpp==1)
                    totalF += world[:tF][q, n] * world[:tF][qOpp, nOpp] * world[:epsilon]
                    # what is focal patch is minimal richness? then can only win losses dont change freq
                    # or if qOpp is maximum and q is less than it
                    if ((q == 1) & (qOpp > 1)) | ((qOpp == world[:q]) & (q < qOpp))
                        # println("only win")
                        partF += world[:tF][q, n] * world[:tF][qOpp, nOpp] * world[:epsilon]
                    # what if focal patch is non-zero but other patch is 0. Then can only lose
                    elseif ((q > 1) & (qOpp == 1)) | ((q == world[:q]) & (qOpp < q))
                        # println("only lose")
                        partF += world[:tF][q, n] * world[:tF][qOpp, nOpp] * world[:epsilon]
                    # if equal richness and neither max not min then they can win or lose 
                    elseif (q == qOpp) & (q > 1) & (q < world[:q])
                        # println("both win or lose")
                        # q wins
                        partF += world[:tF][q, n] * world[:tF][qOpp, nOpp] * world[:epsilon]
                    # or if diffrent and not min or max
                    elseif (q != qOpp) & (q > 1) & (q < world[:q]) & (qOpp > 1) & (qOpp < world[:q])
                        partF += world[:tF][q, n] * world[:tF][qOpp, nOpp] * world[:epsilon]
                    else
                        # println("NO MATCH") 
                    end
                end
                end
            end
        end
    end
    fights = partF/totalF
    return fights
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
        funF = ModelingToolkit.build_function(
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

# w1 = produceSim(worldSet[1])

# world[:nGens] = world[:realGen]
# worldSet = dict_list(world)
ls =[]
for cosm in worldSet
    cosm[:nGens] = cosm[:realGen]
    cosm[:gain] = cosm[:ratio]/cosm[:stab]
    cosm[:loss] = (1-cosm[:ratio])/cosm[:stab]
    @time resWorld = produceSim(cosm)
    finalWorld = testRatios!(resWorld)
    global finalWorld[:fights] = calcFights(finalWorld)
    push!(ls, finalWorld)
    # save(joinpath("..", "data", savename(world, "bson")), resWorld)
end


# using StatsBase
# using Plots 
# gr()

# nArray = repeat([i for i in 1:(world[:n]-1)]', world[:q])
# baseF = []
# baseX = []
# baseY = []
# baseM = []
# rbase = []
# for ww in ls 
#     tF = ww[:tF]
#     iX = mean((tF[:, 2:end] .* ww[:tX][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#     iY = mean((tF[:, 2:end] .* ww[:tY][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#     iM = mean((tF[:, 2:end] .* ww[:M][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#     push!(baseF, sum(tF[:, 2:end].*nArray))
#     push!(baseX, iX)
#     push!(baseY, iY)
#     push!(baseM, iM)
#     push!(rbase, ww[:ratio])
# end


# plt = plot(rbase, baseF.-baseF, label="evolved", title="Population Size", legend= :outerbottom, size=(1800, 1000));
# for i in 1:length(ls)
#     ww = ls[i]
#     Flist = [] 
#     Xlist = []
#     Ylist = []
#     for r in world[:ratio]
#         tF = ww[Symbol(replace(string("rF","_", r), "."=>"_"))]
#         iX = mean(tF[:, 2:end] .* ww[:tX][:, 2:end])
#         iY = mean(tF[:, 2:end] .* ww[:tY][:, 2:end])
#         push!(Flist, sum(tF[:, 2:end].*nArray))
#         push!(Xlist, iX)
#         push!(Ylist, iY)
#     end
#     plt = plot!(plt, world[:ratio], Flist.-baseF, label=string("peturbed ", rbase[i]));
#     plt = scatter!(plt, [rbase[i]], [0], label="");
# end
# @show plt
# savefig(plt, "../graphs/populationSize_peturbed.pdf")

# plt = plot(rbase, baseX.-baseX, label="evolved", title="Cooperation Level", legend=:outerbottom, size=(1800, 1000));
# for i in 1:length(ls)
#     ww = ls[i]
#     Flist = [] 
#     Xlist = []
#     Ylist = []
#     for r in world[:ratio]
#         tF = ww[Symbol(replace(string("rF","_", r), "."=>"_"))]
#         iX = mean((tF[:, 2:end] .* ww[:tX][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#         iY = mean(tF[:, 2:end] .* ww[:tY][:, 2:end])
#         push!(Flist, sum(tF[:, 2:end].*nArray))
#         push!(Xlist, iX)
#         push!(Ylist, iY)
#     end
#     plt = plot!(plt, world[:ratio], Xlist.-baseX, label=string("peturbed ", rbase[i]));
#     plt = scatter!(plt, [rbase[i]], [0], label="");
# end
# @show plt
# savefig(plt, "../graphs/cooperation_peturbed.png")

# plt = plot(rbase, baseY.-baseY, label="evolved", title="Conflict Level", legend=:outerbottom, size=(1800, 1000));
# for i in 1:length(ls)
#     ww = ls[i]
#     Flist = [] 
#     Xlist = []
#     Ylist = []
#     for r in world[:ratio]
#         tF = ww[Symbol(replace(string("rF","_", r), "."=>"_"))]
#         iX = mean(tF[:, 2:end] .* ww[:tX][:, 2:end])
#         iY = mean((tF[:, 2:end] .* ww[:tY][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#         push!(Flist, sum(tF[:, 2:end].*nArray))
#         push!(Xlist, iX)
#         push!(Ylist, iY)
#     end
#     plt = plot!(plt, world[:ratio], Ylist.-baseY, label=string("peturbed ", rbase[i]));
#     plt = scatter!(plt, [rbase[i]], [0], label="");
# end
# @show plt
# savefig(plt, "../graphs/conflict_peturbed.png")

# plt = plot(rbase, baseM.-baseM, label="evolved", title="Mortality Level", legend=:outerbottom, size=(1800, 1000));
# for i in 1:length(ls)
#     ww = ls[i]
#     Flist = [] 
#     Xlist = []
#     Ylist = []
#     Mlist = []
#     for r in world[:ratio]
#         tF = ww[Symbol(replace(string("rF","_", r), "."=>"_"))]
#         iX = mean(tF[:, 2:end] .* ww[:tX][:, 2:end])
#         iY = mean(tF[:, 2:end] .* ww[:tY][:, 2:end])
#         iM = mean((tF[:, 2:end] .* ww[:M][:, 2:end] .* nArray) / (tF[:, 2:end] .* nArray))
#         push!(Flist, sum(tF[:, 2:end].*nArray))
#         push!(Xlist, iX)
#         push!(Ylist, iY)
#         push!(Mlist, iM)
#     end
#     plt = plot!(plt, world[:ratio], Mlist.-baseM, label=string("peturbed ", rbase[i]));
#     plt = scatter!(plt, [rbase[i]], [0], label="");
# end
# @show plt
# savefig(plt, "../graphs/mortality_peturbed.png")

# gradX = substitute.(grads[1], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))
# save("worldList.bson", ls)

# using StatsBase
# testDat = ls[1]
# # testDat = load("tempDicts/1.bson")
# testDat[:relW] = testDat[:tW]./mean(testDat[:tW])
# testDat[:qVal] = mean(mapslices(diff, testDat[:relW], dims=1))
# testDat[:weightedFit] = testDat[:relW] .* testDat[:tF]
# testDat[:weightedQ] =  mean(mapslices(diff, testDat[:weightedFit], dims=1))

# # calc diffs in fitness for each term
# world = copy(testDat)
# totDiff = world[:tW][1,2] - world[:tW][2,2]
# Wn = makeFArray(world)
# Wd = makeFArray(world)

# world = ls[1]
# world[:gain] = world[:ratio]/world[:stab]
# world[:loss] = (1-world[:ratio])/world[:stab]

# world[:tX] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
# world[:tY] = Symbolics.value.(fill!(zeros(world[:q], world[:n]), 1E-6))
# world[:tTr] = simpleTrans(world)

# u = ones(world[:q], world[:n])
# u[:, 1] .= 0.0
# # world[:tW] = u
# world[:tF] = reshape(
#     repeat([1/world[:size]], world[:size]), (world[:q], world[:n])
# )
# # world[:tR] = reshape(
# #     repeat(
# #         hcat([0 0], [1/p for p in 2:(world[:n]-1)]'), 
# #         world[:q]
# #     ), 
# #     (world[:q], world[:n])
# # )
# modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, modelM, modelP, modelC = makeModelExpr(world)

# world = updateMPC(world, modelM, modelP, modelC, modelMl, modelPl, modelCl)

# wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)

# ## grads made in parts 
# grads, directSel, indirectSel = makeSelGrads(Wn./Wd, modelM, modelP, modelC, modelMf, modelMl, modelPf, modelPl, modelCf, modelCl, world);
# # grads = directSel

# gradX = substitute.(grads[1], matSub(
#     X.=>world[:tX], 
#     Y.=>world[:tY], 
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     ))
# gradY = substitute.(grads[2], matSub(
#     X.=>world[:tX], 
#     Y.=>world[:tY], 
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     ))

# gradXdir = substitute.(directSel[1], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))
# gradYdir = substitute.(directSel[2], matSub(X.=>world[:tX], Y.=>world[:tY], F.=>world[:tF], W.=>world[:tW], R.=>world[:tR]))

# ##print a grads 
# using Latexify 
# Xlatex = latexify(grads[1])
# Ylatex = latexify(grads[2])
# numX = latexify(gradX)
# numY = latexify(gradY)
# ##

# ## Grads from X and Y directly 
# # sub in vals to eliminate MPC 
# wSysXY = substitute.(Wn./Wd, matSub(
#     Mf.=>modelMf, 
#     Ml.=>modelMl,
#     P.=>modelP, 
#     Pf.=>modelPf, 
#     Pl.=>modelPl, 
#     C.=>modelC,
#     Cf.=>modelCf, 
#     Cl.=>modelCl, 
#     ))
# wDiffXf = Symbolics.derivative.(wSysXY, Xf)
# wDiffXl = Symbolics.derivative.(wSysXY, Xl)
# wDiffYf = Symbolics.derivative.(wSysXY, Yf)
# wDiffYl = Symbolics.derivative.(wSysXY, Yl)
# inclusiveX = wDiffXf .+ R .* wDiffXl
# inclusiveY = wDiffYf .+ R .* wDiffYl

# X2latex = latexify(inclusiveX);
# Y2latex = latexify(inclusiveY);

# gradX2 = substitute.(
#     substitute.(inclusiveX, 
#         ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],)
#     ), 
#     matSub(
#     Xf.=>world[:tX], 
#     Xl.=>world[:tX],
#     Y.=>world[:tY],
#     Yf.=>world[:tY], 
#     Yl.=>world[:tY],
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     )
# )
# gradY2 = substitute.(
#     substitute.(inclusiveY, 
#         ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],)
#     ), 
#     matSub(
#     Xf.=>world[:tX], 
#     Xl.=>world[:tX],
#     Y.=>world[:tY],
#     Yf.=>world[:tY], 
#     Yl.=>world[:tY],
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     )
# )

# gradX2dir = substitute.(
#     substitute.(wDiffXf, 
#         ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],)
#     ), 
#     matSub(
#     Xf.=>world[:tX], 
#     Xl.=>world[:tX],
#     Y.=>world[:tY],
#     Yf.=>world[:tY], 
#     Yl.=>world[:tY],
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     )
# )
# gradY2dir = substitute.(
#     substitute.(wDiffYf, 
#         ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],)
#     ), 
#     matSub(
#     Xf.=>world[:tX], 
#     Xl.=>world[:tX],
#     Y.=>world[:tY],
#     Yf.=>world[:tY], 
#     Yl.=>world[:tY],
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     )
# )

# gradY2indir = substitute.(
#     substitute.(world[:tR] .* wDiffYl, 
#         ([B=>world[:basem], epsilon=>world[:epsilon], d=>world[:d]],)
#     ), 
#     matSub(
#     Xf.=>world[:tX], 
#     Xl.=>world[:tX],
#     Y.=>world[:tY],
#     Yf.=>world[:tY], 
#     Yl.=>world[:tY],
#     F.=>world[:tF], 
#     W.=>world[:tW], 
#     R.=>world[:tR]
#     )
# )

# wSys[1,2] = 1 -sum((W.*F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))')./sum(F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))'))
# for q in 1:world[:q]
#     wSys[q, 1] = W[q,1]
# end
# wSys12 = copy(wSys)

# wSys, Wn, Wd = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)
# wSys[2,2] = 1 -sum((W.*F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))')./sum(F.*reshape(repeat([x-1 for x in 1:world[:n]], world[:q]), (world[:n],world[:q]))'))
# for q in 1:world[:q]
#     wSys[q, 1] = W[q,1]
# end
# wSys22 = copy(wSys)

# funW = ModelingToolkit.build_function(
#         wSys12, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
#         expression=Val{false}
#         );
# callFunW12 = eval(funW[2]);

# funW = ModelingToolkit.build_function(
#     wSys22, W, F, Mf, Ml, Pf, Pl, P, C, Cl, Cf, d, epsilon;
#     expression=Val{false}
#     );
# callFunW22 = eval(funW[2]);

# fSys = makeFsys(F, M, P, C, d, world)
#     funF = ModelingToolkit.build_function(
#         fSys, F, M, P, C, d, epsilon;
#         expression=Val{false}
#         );
#     callFun = eval(funF[2]);

# wFun12(Fx, x) = callFunW12(
#     reshape(Fx, (world[:q], world[:n])), 
#     reshape(x, (world[:q], world[:n])), 
#     world[:tF], 
#     world[:Mf], 
#     world[:Ml], 
#     world[:Pf], 
#     world[:Pl],
#     world[:P],
#     world[:C],
#     world[:Cl],
#     world[:Cf],
#     world[:d], 
#     world[:epsilon]
# )

# wFun22(Fx, x) = callFunW22(
#     reshape(Fx, (world[:q], world[:n])), 
#     reshape(x, (world[:q], world[:n])), 
#     world[:tF], 
#     world[:Mf], 
#     world[:Ml], 
#     world[:Pf], 
#     world[:Pl],
#     world[:P],
#     world[:C],
#     world[:Cl],
#     world[:Cf],
#     world[:d], 
#     world[:epsilon]
# )

# fFun(Fx, x) = callFun(
#         reshape(Fx, (world[:q], world[:n])), 
#         reshape(x, (world[:q], world[:n])),  
#         world[:M], 
#         world[:P], 
#         world[:C], 
#         world[:d], 
#         world[:epsilon]
#         )
# solF = genSolF(fFun, world)
# world[:tF] = solF.zero


# solW12 = genSolW(wFun12, world)
# # display(solW12.zero)
# solW22 = genSolW(wFun22, world)
# # display(solW22.zero)

# println("Freq")
# display(pprint(world[:tF]))
# println("Fitnesses")
# display(pprint(world[:tW]))
# println("Relatednesses")
# display(pprint(world[:tR]))
# println("X grad from W")
# display(gradX2)
# println("Y grad from W")
# display(gradY2)
# println("X grad from MPC")
# display(gradX)
# println("Y grad from MPC")
# display(gradY)

#jyulia
# -(0.1W₁ˏ₃ + 0.1W₂ˏ₂ + 0.1F₁ˏ₁*W₁ˏ₂ + 0.05F₁ˏ₂*(W₁ˏ₂ + W₁ˏ₃) + 0.05F₂ˏ₁*(W₁ˏ₂ + W₂ˏ₂) + 0.05F₂ˏ₂*(W₁ˏ₂ + W₂ˏ₃) + 0.05F₃ˏ₁*(W₁ˏ₂ + W₃ˏ₂) + 0.05F₃ˏ₂*(W₁ˏ₂ + W₃ˏ₃) + 0.5W₁ˏ₃*(0.1F₁ˏ₂ + 0.2F₁ˏ₃ + 0.4F₂ˏ₂ + 0.5F₂ˏ₃ + 0.7F₃ˏ₂ + 0.8F₃ˏ₃))*((0.15000000000000002 + 0.05F₁ˏ₁ + 0.05F₁ˏ₂ + 0.05F₂ˏ₁ + 0.05F₂ˏ₂ + 0.05F₃ˏ₁ + 0.05F₃ˏ₂ + 0.05F₁ˏ₂ + 0.1F₁ˏ₃ + 0.2F₂ˏ₂ + 0.25F₂ˏ₃ + 0.35F₃ˏ₂ + 0.4F₃ˏ₃ + 0.1exp(-Xf₁ˏ₂) + 0.1(Xf₁ˏ₂^2) + 0.1(Yf₁ˏ₂^2))^-2)*(0.2Xf₁ˏ₂ - (0.1exp(-Xf₁ˏ₂)))

#mathemtica
# -(((0.1* W[1,3]+0.5 *(0.1 *F[1,1]+0.2 *F[1,2]+0.4 *F[2,1]+0.5 *F[2,2]+0.7 *F[3,1]+0.8 *F[3,2]) *W[1,2+1]+0.1 *W[2,1+1]+0.05 *(W[1,1+1]+F[1,0] *W[1,1+1]+F[1,1] *W[1,2+1]+F[2,0] W[2,1+1]+F[2,1] W[2,2+1]+F[3,0] W[3,1+1]+F[3,1] W[3,2+1]))*(-0.1 exp(-Xf[1,1])+0.2 Xf[1,1]))/(0.2 +0.1 exp(-Xf[1,1])+0.5 (0.1 F[1,1]+0.2 F[1,2]+0.4 F[2,1]+0.5 F[2,2]+0.7 F[3,1]+0.8 F[3,2])+0.1 Xf[1,1]^2+0.1 Yf[1,1]^2)^2)

# -(((0.1 w12+0.5 (0.1 f11+0.2 f12+0.4 f21+0.5 f22+0.7 f31+0.8 f32) w12+0.1 w21+0.05 (w11+f10 w11+f11 w12+f20 w21+f21 w22+f30 w31+f31 w32)) (-0.1 E^-xi11+0.2 xi11))/(0.2 +0.1 E^-xi11+0.5 (0.1 f11+0.2 f12+0.4 f21+0.5 f22+0.7 f31+0.8 f32)+0.1 xi11^2+0.1 yi11^2)^2)

# mortality 
# wSys, num, den = makeWsys(W, F, Mf, Ml, P, Pf, Pl, C, Cf, Cl, d, epsilon, world)
# wSys[1,2] = 1 - W[1,2]

# function gradW(wSys, tW, world)
#     return res = substitute.(
#         substitute.(
#             substitute.(
#                 wSys, d=>world[:d]
#             ), 
#             epsilon=>world[:epsilon]
#         ), 
#         matSub(
#             F.=>world[:tF], 
#             Cl.=>world[:Cl], 
#             Cf.=>world[:Cf],
#             C.=>world[:C],
#             Ml.=>world[:Ml],
#             Mf.=>world[:Mf],
#             M.=>world[:M],
#             Pl.=>world[:Pl],
#             Pf.=>world[:Pf],
#             P.=>world[:P],
#             W.=>tW
#         )
#     )
# end 

# Wn1, Wd1 = addMortality(copy(Wn), copy(Wd), W, Mf, Ml, world)

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

# Wn = makeFArray(world)
# Wd = makeFArray(world)
# Wn1, Wd1 = addMortality(copy(Wn), copy(Wd), W, Mf, Ml, world)
# Wn2, Wd2 = addTrans(copy(Wn1), copy(Wd1), W, world)
# Wn3, Wd3 = addLocBirth(copy(Wn2), copy(Wd2), W, Pf, Pl, d, world)
# Wn4, Wd4 = addImmigration(copy(Wn3), copy(Wd3), W, F, P, d, world)
# Wn5, Wd5 = addFights(copy(Wn4), copy(Wd4), W, F, C, Cf, Cl, epsilon, world)
# Wn6, Wd6 = addDistantBirth(copy(Wn5), copy(Wd5), W, F, Pf, world)

# ww = Wn5-Wn4
# ww2 = substitute.(ww, matSub(Cl.=>modelCf, Cf.=>modelCl))
# ww3 = substitute.(ww2, ([epsilon=>world[:epsilon]]))
# ww4 = substitute.(ww3, matSub(Yl.=>world[:tY], Yf.=>world[:tY]))

# wwd = Wd5-Wd4
# wwd2 = substitute.(wwd, matSub(Cl.=>modelCf, Cf.=>modelCl))
# wwd3 = substitute.(wwd2, ([epsilon=>world[:epsilon]]))
# wwd4 = substitute.(wwd3, matSub(Yl.=>world[:tY], Yf.=>world[:tY]))

# ws = Wn6 
# ws1 = substitute.(ws, matSub(
#     [epsilon].=>[1], 
#     [d].=>[0.5],
#     Cl.=>world[:tY],
#     Cf.=>world[:tY],
#     Pf.=>world[:Pf],
#     Ml.=>world[:Ml],
#     Mf.=>world[:Mf]
# ))
# ws1[1,4]

# dws1 = Symbolics.derivative(ws1[2,3], Yl[2,3])

# dw = wDiffYl
# dw = substitute.(dw, matSub(
#     # [epsilon].=>[1], 
#     # [d].=>[0.5],
#     Cl.=>world[:tY],
#     Cf.=>world[:tY],
#     Pf.=>world[:Pf],
#     Ml.=>world[:Ml],
#     Mf.=>world[:Mf]
# ))
# dw[2,3]

# wSysXY = substitute.(Wn6, matSub(
#     Mf.=>modelMf, 
#     Ml.=>modelMl,
#     P.=>modelP, 
#     Pf.=>modelPf, 
#     Pl.=>modelPl, 
#     Cf.=>modelCf, 
#     Cl.=>modelCl, 
#     ))

# dws1 = Symbolics.derivative(wSysXY[2,3], Yl[2,3])

# dws = substitute.(dws1, matSub(
#     [epsilon].=>[1], 
#     [d].=>[0.5],
#     Cl.=>world[:tY],
#     Cf.=>world[:tY],
#     Pf.=>world[:Pf],
#     Ml.=>world[:Ml],
#     Mf.=>world[:Mf]
# ))[1]

# dws = substitute.(Wn6, matSub(
#     [epsilon].=>[1], 
#     [d].=>[0.5],
#     Cl.=>modelCl,
#     Cf.=>modelCf,
#     Pf.=>world[:Pf],
#     Ml.=>world[:Ml],
#     Mf.=>world[:Mf]
# ))[1]

# dws[2,3]

# v = makeFArray(world[:q], world[:n])
# for q in 1:world[:q]
#     for n in 1:world[:n]
#         v[q,n] = victory(Yf[2,3]+ Yl[2,3], (n-1)*Yl[q, n])
#     end
# end

