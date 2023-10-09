using PlotlyJS
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]

# filter df for q=3 n=3
dfplot = filter(df -> df[:q]==3 && df[:n]==3, df)