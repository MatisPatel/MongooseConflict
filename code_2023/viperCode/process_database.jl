using DrWatson
using DataFrames 
using CSVFiles

dat = (joinpath("..", "results_data.jld2"))

df = dat["df"]

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

# make column of pasted togehter mult X and Y 
df[!, :multXY] = [string(x, ":", y) for (x, y) in zip(df[!, :multX], df[!, :multY])]

# make column of pasted togehter shape of X and Y 
df[!, :shapeXY] = [string(x, ":", y) for (x, y) in zip(df[!, :shape_X_cost], df[!, :shape_Y_cost])]

sort!(df, [:ratio, :b, :k, :d])

# save dataframe as jld2 
save("processed_df.jld2", "df", df)