using Symbolics 

for i in 2:10 
    @variables begin 
        X[1:i, 1:i]
        Y[1:i, 1:i]
    end 

    expr = X*Y - X.^2
    expr = Symbolics.scalarize(expr)
    println(expr)
end