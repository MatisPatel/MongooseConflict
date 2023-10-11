using PlotlyJS
using CSVFiles
using DataFrames

df = load("processed_df.jld2")["df"]

# filter df for q=3 n=3
dfplot = filter(df -> df[:q]==3 && 
df[:n]==3 && df[:stab]==2 && df[:shape_X_cost]==2 && df[:shape_Y_cost]==2, df)

plot(dfplot, x=:ratio, y=:meanExprX, color=:epsilon)

plot(dfplot, x=:ratio, y=:meanExprY, color=:epsilon)