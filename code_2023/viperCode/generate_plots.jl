using PlotlyJS
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]
save("plotData.csv", df)

function make_mean_plots(df, q, n)
    dfplot = filter(row -> row[:q] == q && row[:n] == n, df)

    # plot of meanX expression vs ratio and save the plot
    plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanXexpr.png")

    # plot of meanY expression vs harshness and save the plot 
    plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanYexpr.png")

    # plot of meanR vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanR.png")

    # plot of meanW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanW.png")

    # plot of varW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_varW.png")
end

make_mean_plots(df, 3, 3)

# filter df for q=3 n=3
dfplot = filter(df -> 
    # all(df[:tR].>=-0.01) &&
    df[:q]==3 && 
    df[:n]==3 &&
    df[:stab] in [1.0, 3.0] &&
    df[:epsilon] in [1.0, 2.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/3x3_ratio_meanXexpr.png")

# plot of meanY expression vs harshness and save the plot 
plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/3x3_ratio_meanYexpr.png")

# plot of meanR vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/3x3_ratio_meanR.png")

# plot of meanW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/3x3_ratio_meanW.png")

# plot of varW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/3x3_ratio_varW.png")

# filter df for q=4 n=4
dfplot = filter(df -> 
    all(df[:tR].>=-0.01) &&
    df[:q]==4 && 
    df[:n]==4 &&
    df[:stab] in [1.0, 3.0] &&
    df[:epsilon] in [1.0, 2.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x4_ratio_meanXexpr.png")

# plot of meanY expression vs harshness and save the plot 
plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x4_ratio_meanYexpr.png")

# plot of meanR vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x4_ratio_meanR.png")

# plot of meanW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x4_ratio_meanW.png")

# plot of varW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x4_ratio_varW.png")

# filter df for q=5 n=5
dfplot = filter(df -> 
    all(df[:tR].>=-0.01) &&
    df[:q]==5 && 
    df[:n]==5 &&
    df[:stab] in [1.0, 3.0] &&
    df[:epsilon] in [1.0, 2.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/5x5_ratio_meanXexpr.png")

# plot of meanY expression vs harshness and save the plot 
plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/5x5_ratio_meanYexpr.png")

# plot of meanR vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/5x5_ratio_meanR.png")

# plot of meanW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/5x5_ratio_meanW.png")

# plot of varW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/5x5_ratio_varW.png")

# filter df for q=6 n=4
dfplot = filter(df -> 
    all(df[:tR].>=-0.01) &&
    df[:q]==6 && 
    df[:n]==4 &&
    df[:stab] in [1.0, 3.0] &&
    df[:epsilon] in [1.0, 2.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/6x4_ratio_meanXexpr.png")

# plot of meanY expression vs harshness and save the plot 
plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/6x4_ratio_meanYexpr.png")

# plot of meanR vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/6x4_ratio_meanR.png")

# plot of meanW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/6x4_ratio_meanW.png")

# plot of varW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/6x4_ratio_varW.png")

# filter df for q=4 n=6
dfplot = filter(df -> 
    all(df[:tR].>=-0.01) &&
    df[:q]==4 && 
    df[:n]==6 &&
    df[:stab] in [1.0, 3.0] &&
    df[:epsilon] in [1.0, 2.0, 3.0],
df)

# plot of meanX expression vs ratio and save the plot
plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x6_ratio_meanXexpr.png")

# plot of meanY expression vs harshness and save the plot 
plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x6_ratio_meanYexpr.png")

# plot of meanR vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanR, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x6_ratio_meanR.png")

# plot of meanW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:meanW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x6_ratio_meanW.png")

# plot of varW vs harshness and save plot 
plt = plot(dfplot, x=:harshness, y=:varW, color=:epsilon, facet_row=:stab);
savefig(plt, "../graphs/4x6_ratio_varW.png")

