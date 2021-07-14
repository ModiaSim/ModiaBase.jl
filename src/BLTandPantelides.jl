"""
Module with graph theoretical methods for maximum matching, Pantelides index reduction and Tarjan's strongly connected components.

* Author: Hilding Elmqvist, Mogram AB  
* Date: July-August 2016 (rewritten). 
* License: MIT

A bipartite graph (E, V) is defined by an array of arrays. Each E-entry has an integer array of indices to V-entries. 

Example bipartite graph:
 
    G = [
      [3, 5],
      [4, 6],  
      [1, 7, 9],
      [2, 8, 9],
      [1, 2]  
    ]
"""
module BLTandPantelides

export matching, pantelides!, BLT, checkAssign

using ..BLTandPantelidesUtilities

"Controls logging"
const log = false

"""
    function augmentPath!(G, i, assign, vColour, eColour, vPassive)
Construction of augmenting path

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function augmentPath!(G, i, assign, vColour, eColour, vPassive)
    # returns pathFound
    # assign: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
    # i: E-node
    # vPassive: set to != 0 has the same effect as deleting V-node and corresponding edges
    # j: V-node

    if log
        println("augmentPath: equation $i")
    end

    pathFound = false
    eColour[i] = true
    
    # If a V-node j exists such that edge (i-j) exists and assign[j] == 0
    for j in G[i]
        if vPassive[j] == 0 && assign[j] == 0
            pathFound = true
            assign[j] = i
            return pathFound
        end
    end
  
    # For every j such that edge (i-j) exists and j is uncoloured
    for j in G[i]
        if vPassive[j] == 0 && !vColour[j]
            vColour[j] = true
            k = assign[j]
            pathFound = augmentPath!(G, k, assign, vColour, eColour, vPassive)
    
            if pathFound 
                assign[j] = i
                return pathFound
            end
        end
    end
    return pathFound
end


function checkAssign(assign, VSizes, VTypes, ESizes, ETypes, equationsInfix, variableNames, A, vPassive=A)
    println("Checking assignment")
    assignmentOK = true
    for j in 1:length(assign)
        if vPassive[j] == 0
            i = assign[j]
            if i > 0 && VSizes[j] != ESizes[i]
                assignmentOK = false
                print("Error: Variable ") 
                printList(variableNames, [j], A, newLine=false)
                println(" (($j)) with size=$(VSizes[j]) is assigned in equation (($i)) with size $(ESizes[i])")
            end
        end
    end

    if assignmentOK
        println("Assignment is OK")
    else
        # error("Assignment not OK")
  end
end


"""
    function matching(G, M, vActive=fill(true, M))
Find maximum matching in bipartite graph

* `G`: bipartite graph
* `M`: number of V-nodes
* `vActive`: set to false has the same effect as deleting V-node and corresponding edges
* `return assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function matching(G, M, vActive=fill(true, M))
    assign::Array{Int,1} = fill(0, M)
    eColour::Array{Bool,1} = fill(false, length(G))
    vColour::Array{Bool,1} = fill(false, M)
    vPassive::Array{Int,1} = [if va; 0 else 1 end for va in vActive]
    for i in 1:length(G)
        fill!(eColour, false)
        fill!(vColour, false)
        pathFound = augmentPath!(G, i, assign, vColour, eColour, vPassive)
    end
    return assign
end

# -------------------------------------------------------

"""
    function pantelides!(G, M, A)
Perform index reduction with Pantelides algorithm.

* `G`: bipartite graph (updated)
* `M`: number of V-nodes
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return B`: B[i] = if E[l] = der(E[i]) then l else 0

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function pantelides!(G, M, A)
    assign::Array{Int,1} = fill(0, M)
    B::Array{Int,1} = fill(0, length(G))
    eColour::Array{Bool,1} = fill(false, length(G))
    vColour::Array{Bool,1} = fill(false, M)
    N = length(G)
    N2 = N
    for k in 1:N2
        pathFound = false
        i = k
        while !pathFound
            # Delete all V-nodes with A[.] != 0 and all their incidence edges from the graph
            # Designate all nodes as "uncoloured"
            if length(eColour) == length(G)
                fill!(eColour, false)
            else
                eColour = fill(false, length(G))
            end
            if length(vColour) == length(M)
                fill!(vColour, false)
            else
                vColour = fill(false, M)
            end

            pathFound = augmentPath!(G, i, assign, vColour, eColour, A)
            if !pathFound
                if log
                    println("\nDifferentiate:")
                end
                
                # For every coloured V-node j do
                for j in 1:length(vColour)
                    if vColour[j]
                        M += 1
                        if log
                            println("New variable derivative: var($M) = der($j)")
                        end
                        push!(A, 0)
                        A[j] = M
                        push!(assign, 0)
                    end
                end
        
                # For every coloured E-node l do
                for l in 1:N
                    if eColour[l]
                        N += 1
                        if log
                            println("New equation derivative:  equ($N) = DER($l)")
                        end
                        
                        # Create new E-node N
                        push!(G, copy(G[l]))
            
                        # Create edges from E-node N to all V-nodes j and A[j] such that edge (l-j) exists
                        for m in 1:length(G[l])
                            j = G[l][m]
                            if !(A[j] in G[N])
                                push!(G[N], A[j])
                            end
                        end
                        push!(B, 0)
                        
                        # Set B[l] = N
                        B[l] = N
                    end
                end

                # For every coloured V-node j
                for j in 1:length(vColour)
                    if vColour[j]
                        if log
                            println("Assigning derivative of variable $(A[j]) to derivative of equation: $(B[assign[j]])")
                        end
                        assign[A[j]] = B[assign[j]]
                    end
                end
                i = B[i]
            end
        end
    end

    return assign, A, B
end


const notOnStack = 1000000000

"""
Find minimal systems of equations that have to be solved simultaneously.
Reference: 
Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms", SIAM Journal on Computing 1 (2): 146–160, doi:10.1137/0201010
"""
function strongConnect!(G, assign, v, nextnode, stack, components, lowlink, number)
    # println("strongConnect: ", v)
  
    if v == 0 
        return nextnode
    end
    
    nextnode += 1
    lowlink[v] = number[v] = nextnode
    push!(stack, v)

    for w in [assign[j] for j in G[v]] # for w in the adjacency list of v
        if w > 0   # Is assigned
            if number[w] == 0 # if not yet numbered
                nextnode = strongConnect!(G, assign, w, nextnode, stack, components, lowlink, number)
                lowlink[v] = min(lowlink[v], lowlink[w])
            else
                if number[w] < number[v]
                    # (v, w) is a frond or cross-link
                    # if w is on the stack of points. Always valid since otherwise number[w]=notOnStack (a big number)
                    lowlink[v] = min(lowlink[v], number[w])
                end
            end
        end
    end
  
    if lowlink[v] == number[v]
        # v is the root of a component
        # start a new strongly connected component
        comp = []
        repeat = true
        while repeat
            # delete w from point stack and put w in the current component
            # println("delete w from point stack and put w in the current component")
            w = pop!(stack)     
            number[w] = notOnStack
            push!(comp, w)
            repeat = w != v
        end 
        push!(components, comp)
    end
    return nextnode
end


"""
    function BLT(G, assign)

Find Block Lower Triangular structure for a bipartite graph `G` with assignment `assign`
    
* `G`: bipartite graph
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return components`: cell array of components. Each component is a list of indices to E-nodes
"""
function BLT(G, assign)
    nextnode::Int = 0
    stack = []
    components = []
    lowlink = fill(0, length(G))
    number = fill(0, length(G))
  
    for v in 1:length(G)
        if number[v] == 0
            nextnode = strongConnect!(G, assign, v, nextnode, stack, components, lowlink, number)
        end    
    end
    return components
end

end