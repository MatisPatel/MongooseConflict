using DrWatson
using BSON
using DataFrames

datdir = joinpath("..", "data")

files = readdir(datdir)

for file in files 
    try
        row = load(joinpath(datdir, file))
    catch 
        
end

df = collect_results!(datdir)