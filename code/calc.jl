function gain(s, h)
    return(h/(h*s + s))
end

function loss(s, h)
    return(1/(h*s + s))
end

function gl(s, h)
    return([gain(s,h), loss(s, h)])  
end

s = 2:2:20
h = vcat(0.2:0.2:1, 2:2:10)

ls = []
for s1 in s
    for h1 in h 
        push!(ls, gl(s1, h1))
    end
end