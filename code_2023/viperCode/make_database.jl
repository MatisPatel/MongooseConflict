using DrWatson
using BSON
using DataFrames

datdir = joinpath("..", "data")

files = readdir(datdir)

# loop over files and read if any read fails delete the file
for file in files
    try
        load(joinpath(datdir, file))
    catch
        rm(joinpath(datdir, file))
    end
end

df = collect_results!(datdir)