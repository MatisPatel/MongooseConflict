using Random
using StatsBase

mutable struct Group 
    ID::Int
    members::Vector{Ind}
    resources::Float64
end

mutable struct Ind 
    ID::Int
    groupID::Int
    X::Float64 
    Y::Float64
end

function updateResources(groups)
    for group in groups
        group.resources *= 0.9 
        group.resources += 1
    end 
    return groups
end

# define pop params
numGroups = 5
numInds = 1
global maxIndID = 1 

mu = 0.01 
baseReprod = 0.01

# 
groups = [Group(i, [], 1) for i in 1:numGroups];
for group in groups 
    for ind in 1:numInds 
        newInd = Ind(maxIndID, group.ID, rand(), rand())
        maxIndID += 1
        push!(group.members, newInd)
    end
end

