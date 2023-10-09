using DrWatson
using BSON
using DataFrames

datdir = joinpath("..", "data")

df = collect_results!(datdir)