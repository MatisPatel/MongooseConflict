using DrWatson
using BSON
using DataFrames

datdir = joinpath("..", "data", "2024_03_11b")

files = readdir(datdir)

# loop over files and read if any read fails delete the file
for file in files
    try
        println(string("attempting to load file: ", file))
        load(joinpath(datdir, file))
    catch
        println("Error reading file, deleting file from dir")
        rm(joinpath(datdir, file))
    end
end

files = readdir(datdir)

df = collect_results!(datdir)