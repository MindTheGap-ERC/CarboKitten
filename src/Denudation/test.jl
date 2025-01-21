mat = 2 .* ones(3,5,10)

function indi(a::Any)
    for i in CartesianIndices(a)
    a[i] *=  2
    end
    return a
end

function comp(a::Any)
    for i in 1:3
    indi(a[i,:,:])
    end
return a
end

comp(mat)

