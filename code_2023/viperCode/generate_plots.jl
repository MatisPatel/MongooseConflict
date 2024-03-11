using PlotlyJS
using CSVFiles
using DataFrames

input = "processed_results_2024_03_11.jld2"
df = load(joinpath("..", "results", input))["df"]
save("plotData.csv", df)

function make_mean_plots2(df, q, n)
    println("filter df for q = $q and n = $n")
    dfplot = filter(row -> 
    row[:q] == q && 
    row[:n] == n && 
    sum(row[:normF])>0 &&
    row[:multX] == 0.1 &&
    row[:multY] == 0.1 
    , df)

    # plot of meanX expression vs ratio and save the plot
    plt = plot(dfplot, x=:harshness, y=:meanExprX, color=:d, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanXexpr.png")

    # plot of meanY expression vs harshness and save the plot 
    plt = plot(dfplot, x=:harshness, y=:meanExprY, color=:d, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanYexpr.png")

    # plot of meanR vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanR, color=:d, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanR.png")

    # plot of meanW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanW, color=:d, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanW.png")

    # plot of varW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:varW, color=:d, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_varW.png")
end

function make_mean_plots(df, q, n)
    dfplot = filter(row -> 
    row[:q] == q && 
    row[:n] == n && 
    sum(row[:normF])>0
    , df)

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

# group df by q and n 
dfgrouped = groupby(df, [:q, :n])
# make plots for each group
for group in dfgrouped
    q = group[!, :q][1]
    n = group[!, :n][1]
    try
        make_mean_plots2(df, q, n)
    catch
        println("error for q = $q and n = $n")
    end
end




