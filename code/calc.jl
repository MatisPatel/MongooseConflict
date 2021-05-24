function gain(s, h)
    return(h/s)
end

function loss(s, h)
    return((1-h)/s)
end

function gl(s, h)
    return([gain(s,h), loss(s, h)])
end

s = 2:2:20
h = 0.1:0.1:0.9

ls = []
for s1 in s
    for h1 in h 
        push!(ls, gl(s1, h1))
    end
end