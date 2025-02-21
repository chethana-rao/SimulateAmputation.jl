using Comodo


function fixteeth(F,B; method=:add)

    if method == :add
        # Add interior teet 
        indNodesOut = unique(reduce(vcat,F[B]))
        B_fixed = Vector{Bool}(undef,length(B))
        for i in eachindex(B)
            B_fixed[i] = all([in(i,indNodesOut) for i in F[i]])
        end
    elseif method == :remove
        indNodesOut = unique(reduce(vcat,F[.!B]))
        B_fixed = Vector{Bool}(undef,length(B))
        for i in eachindex(B)
            B_fixed[i] = !all([in(i,indNodesOut) for i in F[i]])
        end
    else(throw(ArgumentError("Method should be :either add or :remove ")))
    
    end
    return B_fixed
end
