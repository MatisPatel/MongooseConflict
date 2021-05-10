using BSON  
using DataFrames
using DrWatson
using CSV

datdir = joinpath("..", "data")
resdir = joinpath("..", "results")

files  = readdir(datdir)

rows=[]
fullDat = DataFrame(Dict(:ID => 0))
for i in 1:5
    testDat = load(joinpath(datdir, files[i]))
    normF = testDat[:tF][:, 2:end]./sum(testDat[:tF][:, 2:end])
    nArray = repeat([i for i in 1:(testDat[:n]-1)]', testDat[:q])
    testDat[:tXw] = testDat[:tX][:, 2:end] .* normF
    testDat[:tYw] = testDat[:tY][:, 2:end] .* normF
    testDat[:groupAvgX] = sum(testDat[:tXw])
    testDat[:indAvgX] = sum(testDat[:tXw] .*  nArray)/ sum(normF .* nArray)
    testDat[:groupAvgY] = sum(testDat[:tYw])
    testDat[:indAvgY] = sum(testDat[:tYw] .*  nArray)/ sum(normF .* nArray)

    tempDict = Dict{Symbol, Any}(:ID=>i)
    for (key, val) in testDat
        if !isa(val, Array)
            tempDict[key] = val 
        else 
            for q in 1:size(val)[1]
                for n in 1:size(val)[2]
                    tempDict[Symbol(key, q, n-1)] = val[q, n]
                    tempDict[:tn] = n-1 
                    tempDict[:tq] = q
                end
            end
        end
    end
    save(joinpath(resdir, string(tempDict[:ID], ".bson")), tempDict)
    # rowDat = DataFrame(;tempDict...)
    # global fullDat = outerjoin(fullDat, rowDat, on = names(rowDat))
    # append!(fullDat, rowDat)
    # push!(rows, rowDat)
end

# make CSV 
df = collect_results(resdir)
CSV.write(joinpath(resdir, "testDat"), df)