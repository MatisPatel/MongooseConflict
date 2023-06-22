using Symbolics 

@variables begin 
    F[1:2, 1:2]
    W[1:2, 1:2]
end

sys = Array{Num}(undef, 2, 2)

for q in 1:2
    for n in 1:2
        sys[q, n] = F[q,n]*rand() + W[q, n]*rand() 
    end 
end 

Fp = Symbolics.scalarize(W'W)

substitute(sys, Dict([W => Fp]))

