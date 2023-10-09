using BSON
using DrWatson 
using DataFrames
using StatsBase
# using PlotlyJS
# using Tidier

datdir = joinpath("..", "data")

df = collect_results!(datdir)
q = first(df.q)

# add relW col to df 
df[!, :relW] = df[!, :tW]./mean.(df[!, :tW])
# add normalised tF as 1st col are states with no inds 
truncatedF = [x[:, 2:end] for x in df.tF]
normF = [x./sum(x) for x in truncatedF]
df[!, :normF] = normF 

# add truncated tX and tY as 1st col are states with no inds
truncatedX = [x[:, 2:end] for x in df.tX]
truncatedY = [x[:, 2:end] for x in df.tY]
df[!, :tXtrunc] = truncatedX
df[!, :tYtrunc] = truncatedY

# get expressed tX and tY by multiplying truncated X and Y by normF 
expressedX = [x .* y for (x, y) in zip(truncatedX, normF)]
expressedY = [x .* y for (x, y) in zip(truncatedY, normF)]
df[!, :tXexpr] = expressedX
df[!, :tYexpr] = expressedY

#  get overall means of expressed X and Y and variances 
df[!, :meanExprX] = [mean(x) for x in expressedX]
df[!, :meanExprY] = [mean(x) for x in expressedY]
df[!, :varExprX] = [var(x) for x in expressedX]
df[!, :varExprY] = [var(x) for x in expressedY]
df[!, :factor] = [string(x, ":", y, ":", z) for (x, y, z) in zip(df[!, :k], df[!, :b], df[!, :d])]

dfplot = filter(df -> 
df[:epsilon] == 2 && df[:q]==5 && all(df[:tR].>-0.01) && all(df[:tR].<=1), df)
sort!(dfplot, [:ratio, :b, :k, :d])
# plot mean expressed X and Y against ratio
# Plots.plot(
#     dfplot[!, :ratio], 
#     dfplot[!, :meanExprX], 
#     group = dfplot[!, :factor], 
#     label = permutedims(unique(dfplot[!, :factor])),
#     linestyle = :solid, 
#     xlabel = "ratio", 
#     ylabel = "mean expressed X", 
#     title = "mean expressed X against ratio", 
#     legend = :outerbottom,
#     size = (1080, 1900))
# savefig("../results/meanExprX.png")

# Plots.plot(
#     dfplot[!, :ratio], 
#     dfplot[!, :meanExprY], 
#     group = dfplot[!, :factor], 
#     label = permutedims(unique(dfplot[!, :factor])), 
#     linestyle = :solid,
#     xlabel = "ratio", 
#     ylabel = "mean expressed Y", 
#     title = "mean expressed Y against ratio",
#     legend = :outerbottom,
#     size = (1080, 1900))
# savefig("../results/meanExprY.png")

# plt = PlotlyJS.plot(
#     dfplot,
#     x = :ratio, 
#     y = :meanExprY,
#     color = :k,
#     facet_row = :d,
#     facet_col = :b,
# ); 

# plt = PlotlyJS.plot(
#     dfplot,
#     x = :ratio, 
#     y = :meanExprX,
#     color = :k,
#     facet_row = :d,
#     facet_col = :b,
# );



# key_list = [
#     :q,
#     :n,
#     :decay, 
#     :epsilon, 
#     :learning_rate, 
#     :ratio, 
#     :stab]
# for file in readdir(datdir)
#     # load each file and print name 
#     println(file)
#     data = wload(joinpath(datdir, file))
#     shortName = savename(data, "bson", accesses=key_list)
#     println(shortName)
#     wsave(joinpath(datdir, shortName), data)
# end

# plot(dfplot, x=:ratio, y=:meanExprY, color=:d)