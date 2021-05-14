using BSON  
using DataFrames
using DrWatson
using CSV
using StatsBase

datdir = joinpath("..", "data")
resdir = joinpath("..", "results")

files  = readdir(datdir)

function mortFun(n, x, y, B, multX, multY)
    calc = B*2.718^(-1 * (n)*((x*(n-1) + x)/n)) + multX*x^2 + multY*y^2
    return calc
end

rows=[]
fullDat = DataFrame(Dict(:ID => 0))
for i in 1:5
    testDat = load(joinpath(datdir, files[i]))
    testDat[:avgR] = mean(testDat[:tR][:, 3:end])
    normF = testDat[:tF][:, 2:end]./sum(testDat[:tF][:, 2:end])
    nArray = repeat([i for i in 1:(testDat[:n]-1)]', testDat[:q])
    testDat[:tXw] = testDat[:tX][:, 2:end] .* normF
    testDat[:tYw] = testDat[:tY][:, 2:end] .* normF
    testDat[:groupAvgX] = sum(testDat[:tXw])
    testDat[:indAvgX] = sum(testDat[:tXw] .*  nArray)/ sum(normF .* nArray)
    testDat[:groupAvgY] = sum(testDat[:tYw])
    testDat[:indAvgY] = sum(testDat[:tYw] .*  nArray)/ sum(normF .* nArray)

    AnArray = nArray[:, 1]
    testDat[:AtXw] = testDat[:tXw][:, 1:1]
    testDat[:AtYw] = testDat[:tYw][:, 1:1]
    testDat[:AgroupAvgX] = sum(testDat[:AtXw])
    testDat[:AindAvgX] = sum(testDat[:AtXw] .*  AnArray)/ sum(normF .* nArray)
    testDat[:AgroupAvgY] = sum(testDat[:AtYw])
    testDat[:AindAvgY] = sum(testDat[:AtYw] .*  AnArray)/ sum(normF .* nArray)

    SnArray = nArray[:, 2:end]
    testDat[:StXw] = testDat[:tX][:, 3:end] .* normF[:, 2:end]
    testDat[:StYw] = testDat[:tY][:, 3:end] .* normF[:, 2:end]
    testDat[:SgroupAvgX] = sum(testDat[:StXw])
    testDat[:SindAvgX] = sum(testDat[:StXw] .*  SnArray)/ sum(normF .* nArray)
    testDat[:SgroupAvgY] = sum(testDat[:StYw])
    testDat[:SindAvgY] = sum(testDat[:StYw] .*  SnArray)/ sum(normF .* nArray)

    testDat[:avgMort] = mean(mortFun.(nArray, testDat[:tX][:,2:end], testDat[:tY][:,2:end], testDat[:basem], testDat[:multX], testDat[:multY]))
    testDat[:meanFit] = mean(testDat[:tW][:, 2:end])
    testDat[:fit1] = mean(testDat[:tW][1, 2:end])
    testDat[:fit2] = mean(testDat[:tW][2, 2:end])
    testDat[:qVal] = mean(mapslices(diff, testDat[:tW], dims=1))
 
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
CSV.write(joinpath(resdir, "testDat.csv"), df)