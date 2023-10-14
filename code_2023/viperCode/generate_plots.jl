using PlotlyJS
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]

# filter df for q=3 n=3
dfplot = filter(df -> 
    all(df[:tR].>=-0.1) &&
    df[:q]==5 && 
    df[:n]==5 && 
    # df[:epsilon] == 2 &&
    df[:stab]==2 && 
    df[:shape_X_cost]==2 && 
    df[:shape_Y_cost]==2, 
df)

acc = zeros(first(dfplot.q),first(dfplot.n))
n = length(dfplot.tR)
for x in dfplot.tR 
    acc += x./n
end
display(acc)
println(mean(acc[:, 3:end]))

plot(dfplot, x=:ratio, y=:meanExprX, color=:q)

plot(dfplot, x=:ratio, y=:meanExprY, color=:q)

plot(dfplot, x=:ratio, y=:meanR, color=:epsilon)

plot(dfplot, x=:ratio, y=:popSize, color=:epsilon)

# PlotlyJS plot histogram of :timings 
plot(dfplot, x=:timing, color=:epsilon, kind=:histogram)