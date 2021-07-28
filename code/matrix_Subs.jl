using Symbolics

@variables x, y, z 
@variables A[1:3,1:3], B[1:3,1:3], C[1:3,1:3]

scalar = x+y^2 - z 
s1 = substitute(scalar, x=>y)
s2 = substitute(scalar, Dict([x=>y, z=>x]))

matrix = A .+ B.^2 .- C
m1 = substitute(matrix, A=>B)
m2 = substitute(matrix, Dict([A.=>B, C.=>A]))
# 3×3 Matrix{SymbolicUtils.Add{Real, Int64, Dict{Any, Number}, Nothing}}:
#  A[1, 1] + B[1, 1]^2 - C[1, 1]  B[1, 2]^2 + A[1, 2] - C[1, 2]  A[1, 3] + B[1, 3]^2 - C[1, 3]
#  B[2, 1]^2 + A[2, 1] - C[2, 1]  A[2, 2] + B[2, 2]^2 - C[2, 2]  A[2, 3] + B[2, 3]^2 - C[2, 3]
#  A[3, 1] + B[3, 1]^2 - C[3, 1]  B[3, 2]^2 + A[3, 2] - C[3, 2]  B[3, 3]^2 + A[3, 3] - C[3, 3]

function matSub(pairs...)
    pairList = []
    for pair in pairs
        push!(pairList, vec(reshape(pair, (length(pair), 1))))
    end
    return Dict(collect(Iterators.flatten((pairList))))
end

