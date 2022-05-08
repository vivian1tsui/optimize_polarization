include("Graph.jl")
include("GraphMatrix.jl")
include("Algorithm.jl")

function createQ(G :: Graph, n2, Methold) # select n2 point with different methold. || Methold == 'r'(random) | 'l'(largest degree) | 's'(smallest degree)
	if 'r' in Methold
		Q = RandomSeq(G.n, n2)
	else
		D = getDegree(G)
		sD = Array{Tuple{Int32, Float64}}(G.n)
		for i = 1 : G.n
			sD[i] = (i, D[i])
		end
		sort!(sD, by = x -> x[2])
		if 'l' in Methold
			reverse!(sD)
		end
		Q = []
		for i = 1 : n2
			(u, du) = sD[i]
			push!(Q, u)
		end
	end
	f1 = open("Q.txt", "w")
	for x in Q
		println(f1, G.V[x])
	end
	return Q
end

function readQ(G :: Graph)
	d = Dict{Int32, Int32}()
	for i = 1 : G.n
		d[G.V[i]] = i
	end
	Q = []
	open("Q.txt") do f1
		for line in eachline(f1)
			buf = split(line)
			x = parse(Int32, buf[1])
			push!(Q, d[x])
		end
	end
	return Q
end

function createEQ(G :: Graph, Q) # create an edge set EQ while Q is given
	avgw = getAvg(G)
	A = getSparseA(G)
	h = Array{Int8}(G.n)
	fill!(h, 0)
	for x in Q
		h[x] = 1
	end
	VQ = []
	for i = 1 : G.n
		if h[i] == 0
			push!(VQ, i)
		end
	end
	EQ = []
	for x in Q
		for y in VQ
			if A[x, y] == 0
				push!(EQ, (y, x, avgw))
			end
		end
	end
	return EQ
end