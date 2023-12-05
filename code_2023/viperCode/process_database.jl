using DrWatson
using DataFrames 
using CSVFiles
using StatsBase

# function to calculate means and vars for each peterbed run 
function calc_mean_vars(ratio, df)
    return nothing 
end

dat = load(joinpath("..", "data", "results_2023_11_29.jld2"))

df = dat["df"]
df = filter(df-> !(isnan(df[:err])), df)
# add relW col to df 
df[!, :relW] = df[!, :tW]./mean.(df[!, :tW])
# add normalised tF as 1st col are states with no inds 
truncatedF = [x[:, 2:end] for x in df.tF]
normF = [x./sum(x) for x in truncatedF]
df[!, :normF] = normF 
df[!, :popSize] = [sum(x) for x in truncatedF]
df[!, :harshness] = [1-x for x in df.ratio]

# add truncated tW as 1st col are states with no inds
truncatedW = [x[:, 2:end] for x in df.tW]
df[!, :tWtrunc] = truncatedW
# add var and mean of truncated tW
df[!, :varW] = [var(x) for x in truncatedW]
df[!, :meanW] = [mean(x) for x in truncatedW]

# add truncated tX and tY as 1st col are states with no inds
truncatedX = [x[:, 2:end] for x in df.tX]
truncatedY = [x[:, 2:end] for x in df.tY]
df[!, :tXtrunc] = truncatedX
df[!, :tYtrunc] = truncatedY

# get truncated tR as 1st and second columns are not defined 
truncatedR = [x[:, 3:end] for x in df.tR]
df[!, :tRtrunc] = truncatedR
# get mean relatedness 
df[!, :meanR] = [mean(x) for x in truncatedR]

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

# make column of pasted togehter mult X and Y 
df[!, :multXY] = [string(x, ":", y) for (x, y) in zip(df[!, :multX], df[!, :multY])]

# make column of pasted togehter shape of X and Y 
df[!, :shapeXY] = [string(x, ":", y) for (x, y) in zip(df[!, :shape_X_cost], df[!, :shape_Y_cost])]

sort!(df, [:q, :n, :ratio, :b, :k, :d])

# save dataframe as jld2 
save("processed_df.jld2", "df", df)