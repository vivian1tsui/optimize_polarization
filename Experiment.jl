include("Graph.jl")
include("GraphMatrix.jl")
include("Algorithm.jl")
include("SetCreation.jl")

function runExperiment(G :: Graph, Q, EQ, k, alg)
	TIME = time()
	S = []
	if 'O' in alg
		S = Optimal(G, Q, EQ, k)
	elseif 'R' in alg
		S = RandomSelect(EQ, k)
	elseif 'C' in alg
		S = TopCent(G, Q, EQ, k)
	elseif 'D' in alg
		S = TopDegree(G, Q, EQ, k)
	elseif 'E' in alg
		S = ExactSM(G, Q, EQ, k)
	else
		ss = split(alg, ':')
		eps = parse(Float64, ss[2])
		S = ApproxiSM(G, Q, EQ, k, eps)
	end
	TIME = time() - TIME
	return S, TIME
end

function printResult(G :: Graph, Q, io, algs, ALLS, ALLT)
	Total = size(algs, 1)
	for i = 1 : Total
		print(io, algs[i], " ", ALLT[i])
		if (Total == 1) && ('A' in algs[1])
			println(io)
			continue
		end
		if 'O' in algs[i]
			r = EffectiveResistance(G, Q, ALLS[i])
			println(io, " ", r)
			continue
		end
		k = size(ALLS[i], 1)
		tS = []
		for j = 1 : k
			push!(tS, ALLS[i][j])
			r = EffectiveResistance(G, Q, tS)
			print(io, " ", r)
		end
		println(io)
	end
end