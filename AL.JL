include("Graph.jl")
include("GraphMatrix.jl")
using Laplacians

function EffectiveResistance(G :: Graph, Q, S) # Calculate R_q(S)
	L, label = getLQ(G, Q)
	for (u, v, w) in S
		L[label[u], label[u]] += w
	end
	Lp = inv(L)
	return trace(Lp)
end

function combinationNumber(n, k) # Calculate C(n,k)
	ans = 1
	for i = 1 : k
		ans *= (n-i+1)
	end
	for i = 2 : k
		ans /= i
	end
	return ans
end

function nextCombination(A, n, k) # Return next Combination of C(n,k) while current combination is A
	B = A
	for i = k : -1 : 1
		if B[i] < (n-k+i)
			B[i] = B[i] + 1
			for j = i+1 : k
				B[j] = B[j-1]+1
			end
			break
		end
	end
	return B
end

function Optimal(G :: Graph, Q, EQ, k) # Optimal Algorithm
	n = size(EQ, 1)
	tot = combinationNumber(n, k)
	A = Array{Int32}(k)
	for i = 1 : k
		A[i] = i
	end
	opt = 1e10
	optS = []
	for i = 1 : tot
		S = []
		for x in A
			push!(S, EQ[x])
		end
		cnt = EffectiveResistance(G, Q, S)
		if (cnt < opt)
			opt = cnt
			optS = S
		end
		A = nextCombination(A, n, k)
	end
	return optS
end

function RandomSeq(n, k) # Select k numbers randomly from n numbers
	srand(Int(round(time() * 10000.0)))
	A = randcycle(n)
	B = []
	for i = 1 : k
		push!(B, A[i])
	end
	return B
end

function RandomSelect(EQ, k) # Random Algorithm
	n = size(EQ, 1)
	A = randcycle(n)
	S = []
	for i = 1 : k
		push!(S, EQ[A[i]])
	end
	return S
end

function getDegree(G :: Graph) # Return a degree array D
	D = zeros(G.n)
	for (ID, u, v, w) in G.E
		D[u] += w
		D[v] += w
	end
	return D
end

function TopDegree(G :: Graph, Q, EQ, k) # TopDegree Algorithm
	D = getDegree(G)
	for x in Q
		D[x] = -1.0
	end
	tmp = size(EQ, 1)
	h = zeros(tmp)
	S = []
	for i = 1 : k
		md = 0
		select = 0
		for j = 1 : tmp
			(u, v, w) = EQ[j]
			if (h[j] == 0) && (D[u] > md)
				md = D[u]
				select = j
			end
		end
		h[select] = 1
		push!(S, EQ[select])
	end
	return S
end

function TopCent(G :: Graph, Q, EQ, k) # TopCent Algorithm
	Lp = getLp(G)
	for x in Q
		Lp[x, x] = -1e10
	end
	tmp = size(EQ, 1)
	h = zeros(tmp)
	S = []
	for i = 1 : k
		md = -1e9
		select = 0
		for j = 1 : tmp
			(u, v, w) = EQ[j]
			if (h[j] == 0) && (Lp[u, u] > md)
				md = Lp[u, u]
				select = j
			end
		end
		h[select] = 1
		push!(S, EQ[select])
	end
	return S
end

function ExactSM(G :: Graph, Q, EQ, k) # Algorithm 1 in paper
	Lp, label = getLQp(G, Q)
	n = size(Lp, 1)
	S = []
	tot = size(EQ, 1)
	h = zeros(tot)
	fz = zeros(n)
	for i = 1 : k
		# Calculate the delta
		for j = 1 : n
			fz[j] = norm(Lp[:, j])^2
		end
		# Select the edge
		md = 0.0
		select = 0
		for j = 1 : tot
			(u, v, w) = EQ[j]
			u = label[u]
			d = (w * fz[u]) / (1.0 + w * Lp[u, u])
			if (d > md) && (h[j] == 0)
				md = d
				select = j
			end
		end
		h[select] = 1
		push!(S, EQ[select])
		# Update the matrix
		(u, v, w) = EQ[select]
		Lp = updateMatrix(Lp, label[u], w)
	end
	return S
end

function GainsEST(G :: Graph, label, Lsq, EQ, eps) # Algorithm 2 in paper
	srand(Int(round(time() * 10000.0)))
	# Set delta and M
	minw, maxw = getMinMax(G)
	delta1 = (minw^2.0 * eps) / (36.0 * G.n^1.5 * maxw^4.0)
	delta2 = sqrt(eps * minw^2.5 / G.n) / (6.0 * G.n^0.5 * maxw)
	M = round(Int32, 24*log2(2.0*G.n)/(eps^2))
	# Estimate the trace
	f = approxCholSddm(Lsq, tol=1e-6)
	n = size(Lsq, 1)
	ti = zeros(n)
	for i = 1 : M
		z = randn(n, 1)
		y = f(z[:, 1])
		for j = 1 : n
			ti[j] += (y[j, 1]^2)
		end
	end
	for i = 1 : n
		ti[i] /= M
	end
	# Estimate the effective resistance
	##f = approxCholSddm(Lsq, tol=delta2)
	k = round(Int32, 4.0/((eps^2)/2.0+(eps^3)/3.0) * log2(G.n))
	B = getB(G, label)
	X = getX(Lsq)
	m = size(B, 1)
	fm = zeros(n)
	for i = 1 : k
		q = randn(m, 1)
		r = randn(n, 1)
		q1 = B'*q
		r1 = X*r
		z2 = f(q1[:, 1])
		z3 = f(r1[:, 1])
		for j = 1 : n
			fm[j] += (z2[j, 1]^2 + z3[j, 1]^2)
		end
	end
	for i = 1 : n
		fm[i] /= k
	end
	# Calculate the answer
	tot = size(EQ, 1)
	Delta = Array{Float64}(tot)
	for i = 1 : tot
		(u, v, w) = EQ[i]
		Delta[i] = (w*ti[label[u]]) / (1+w*fm[label[u]])
	end
	return Delta
end

function ApproxiSM(G :: Graph, Q, EQ, k, eps) # Algorithm 3 in paper
	S = []
	tot = size(EQ, 1)
	h = zeros(tot)
	sL, label = getSparseLQ(G, Q)
	for i = 1 : k
		Delta = GainsEST(G, label, sL, EQ, eps)
		md = 0
		select = 0
		for j = 1 : tot
			if (Delta[j] > md) && (h[j] == 0)
				md = Delta[j]
				select = j
			end
		end
		h[select] = 1
		push!(S, EQ[select])
		(u, v, w) = EQ[select]
		sL[label[u], label[u]] += w
	end
	return S
end

