using PlotlyJS
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]

plot(df, x=:ratio, y=:meanExprX, color=:shapeXY)