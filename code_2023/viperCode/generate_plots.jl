using Plots
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]

# acc = zeros(first(dfplot.q),first(dfplot.n))
# n = length(dfplot.tR)
# for x in dfplot.tR 
#     acc += x./n
# end
# display(acc)
# println(mean(acc[:, 3:end]))

# filter df for q=3 n=3
dfplot = filter(df -> 
    all(df[:tR].>=-0.01) &&
    df[:q]==3 && 
    df[:n]==3 &&
    df[:stab] in [1.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:ratio, y=:meanExprX, color=:epsilon, facet_row=:stab,
layout=Layout(title = "Mean X expression vs ratio", 
    xaxis_title = "Ratio", 
    yaxis_title = "Mean X expression"))
savefig(plt, "../graphs/3x3_ratio_meanXexpr.png")

plot(dfplot, x=:ratio, y=:meanExprY, color=:epsilon)

plot(dfplot, x=:ratio, y=:meanR, color=:epsilon)

plot(dfplot, x=:ratio, y=:popSize, color=:epsilon)

# PlotlyJS plot histogram of :timings 
plot(dfplot, x=:timing, color=:epsilon, kind=:histogram)