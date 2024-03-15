using PlotlyJS
using CSVFiles
using DataFrames

input = "processed_results_2024_03_11b.jld2"
df = load(joinpath("..", "results", input))["df"]
save("plotData.csv", df)

function make_mean_plots(df, q, n, colorCol=:d)
    println("filter df for q = $q and n = $n")
    dfplot = filter(row -> 
    row[:q] == q && 
    row[:n] == n && 
    sum(row[:normF])>0 &&
    row[:multX] == 0.1 &&
    row[:multY] == 0.1 &&
    row[:peturbed] == false
    , df)

    # plot of meanX expression vs ratio and save the plot
    plt = plot(dfplot, x=:harshness, y=:meanExprX, color=colorCol, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanXexpr.png")

    # plot of meanY expression vs harshness and save the plot 
    plt = plot(dfplot, x=:harshness, y=:meanExprY, color=colorCol, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanYexpr.png")

    # plot of meanR vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanR, color=colorCol, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanR.png")

    # plot of meanW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:meanW, color=colorCol, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_meanW.png")

    # plot of varW vs harshness and save plot 
    plt = plot(dfplot, x=:harshness, y=:varW, color=colorCol, facet_row=:k, facet_col=:b);
    savefig(plt, "../graphs/$(q)x$(n)_ratio_varW.png")
end

# group df by q and n 
dfgrouped = groupby(df, [:q, :n])
# make plots for each group
for group in dfgrouped
    q = group[!, :q][1]
    n = group[!, :n][1]
    try
        make_mean_plots(df, q, n, :basem)
    catch
        println("error for q = $q and n = $n")
    end
end

trace_list = []

for real_ratio in 0.1:0.1:0.9
# try to add plots about peturbations 
dfplot = filter(row -> 
    row[:q] == 4 && 
    row[:n] == 4 && 
    sum(row[:normF])>0 &&
    row[:multX] == 0.1 &&
    row[:multY] == 0.1 &&
    row[:basem] == 0.1 &&
    row[:ratio] == real_ratio
    , df)

tF = Vector{Matrix{Float64}}(undef, 9)
tY = Vector{Matrix{Float64}}(undef, 9)
for (i, ratio) in enumerate(0.1:0.1:0.9)
    tF[i] = dfplot[!, string("tF_", ratio)][i]
    tY[i] = dfplot[!, "tY"][i]
end
truncatedF = [x[:, 2:end] for x in tF]
normF = [x./sum(x) for x in truncatedF]
truncatedY = [x[:, 2:end] for x in tY]
expressedY = [x .* y for (x, y) in zip(truncatedY, normF)]
meanExprY = [mean(x) for x in expressedY]
offset_meanExprY = meanExprY .- dfplot[!, :meanExprY][1:9]
push!(trace_list, offset_meanExprY)
end

matTraces = reduce(hcat, trace_list)

plot(matTraces)