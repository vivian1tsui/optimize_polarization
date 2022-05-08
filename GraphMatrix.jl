include("Graph.jl")

function getL(G :: Graph) # Get Laplacian Matrix
    L :: Matrix{Float32} = zeros(G.n, G.n)
    for (ID, u, v, w) in G.E
        L[u, v] -= w
        L[v, u] -= w
        L[u, u] += w
        L[v, v] += w
    end
    return L
end

function getSparseA(G :: Graph) # Get Sparse Adjacency Matrix
    A :: SparseMatrixCSC{Float64} = spzeros(G.n,G.n)
    for (ID, u, v, w) in G.E
        A[u, v] += w
        A[v, u] += w
    end
    return A
end 

function pseudoInv(A) # Get the Pseudoinverse
    n = size(A, 1)
    Ap = A
    Ap .-= (1.0 / n)
    Ap = inv(Ap)
    Ap .+= (1.0 / n)
    return Ap
end

function getLp(G :: Graph) # Get the Pseudoinverse of Laplacian Matrix
    L = getL(G)
    Lp = pseudoInv(L)
    return Lp
end

function getLQ(G :: Graph, Q) # Calculate L(V \ Q)
    label = Array{Int32}(G.n)
    fill!(label, 0)
    for x in Q
        label[x] = -1
    end
    newSize = 0
    for i = 1 : G.n
        if label[i] == 0
            newSize = newSize + 1
            label[i] = newSize
        end
    end
    L :: Matrix{Float32} = zeros(newSize, newSize)
    for (ID, u, v, w) in G.E
        u1 = label[u]
        v1 = label[v]
        if u1 > 0
            L[u1, u1] += w
        end
        if v1 > 0
            L[v1, v1] += w
        end
        if (u1 > 0) && (v1 > 0)
            L[u1, v1] -= w
            L[v1, u1] -= w
        end
    end
    return L, label
end

function getLQp(G :: Graph, Q) # Calculate L(V \ Q)^-1
    L, label = getLQ(G, Q)
    lu(L)
    svd(L)
    qr(L)
    qr(L)
    qr(L)
    qr(L)
    qr(L)
    Lp = inv(L)
    return Lp, label
end

function getSparseLQ(G :: Graph, Q) # Calculate Sparse Matrix L(V \ Q)
    label = Array{Int32}(G.n)
    fill!(label, 0)
    for x in Q
        label[x] = -1
    end
    newSize = 0
    for i = 1 : G.n
        if label[i] == 0
            newSize = newSize + 1
            label[i] = newSize
        end
    end
    sL :: SparseMatrixCSC{Float64} = spzeros(newSize, newSize)
    for (ID, u, v, w) in G.E
        u1 = label[u]
        v1 = label[v]
        if u1 > 0
            sL[u1, u1] += w
        end
        if v1 > 0
            sL[v1, v1] += w
        end
        if (u1 > 0) && (v1 > 0)
            sL[u1, v1] -= w
            sL[v1, u1] -= w
        end
    end
    return sL, label
end

function getB(G :: Graph, label) # Get Matrix B
    n = 0
    for x in label
        if x > 0
            n = n + 1
        end
    end
    EB = []
    for (ID, u, v, w) in G.E
        if (label[u] > 0) && (label[v] > 0)
            push!(EB, (label[u], label[v], sqrt(w)))
        end
    end
    m = size(EB, 1)
    B = spzeros(m, n)
    for i = 1 : m
        (u, v, w) = EB[i]
        B[i, u] = sqrt(w)
        B[i, v] = -sqrt(w)
    end
    return B
end

function getX(Lsq) # Get Matrix X^0.5
    n = size(Lsq, 1)
    b = zeros(n, 1)
    fill!(b, 1.0)
    bb = Lsq * b
    X = spzeros(n, n)
    for i = 1 : n
        X[i, i] = sqrt(bb[i, 1])
    end
    return X
end

function updateMatrix(A, b, c) # Calculate (B+c*e_b*e_b')^-1 while know B^-1 == A
    bb = A[:, b]
    su = c / (1.0+A[b, b]*c)
    n = size(A, 1)
    for i = 1 : n
        for j = 1 : n
            A[i, j] = A[i, j] - (bb[i, 1] * bb[j, 1] * su) 
        end
    end
    return A
end
