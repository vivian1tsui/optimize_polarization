using Random
using LinearAlgebra
using SparseArrays

function getResult(n, E)
	L = zeros(n, n)
	for (u, v) in E
		L[u, u] += 1.0
		L[v, v] += 1.0
		L[u, v] -= 1.0
		L[v, u] -= 1.0
	end
	Eig = eigvals(L)
	s = 0.0
	for i = 2 : n
		s += (Eig[i]^(-2))
	end
	s /= (2*n)
	return s
end

rng = MersenneTwister(Int(round(time() * 10000.0)))

rlt = spzeros(20000)
for i = 1 : 20
	rlt[i*1000] = 1
end

edge = []
push!(edge, (1, 2))
push!(edge, (1, 3))
push!(edge, (2, 3))

ss = []
push!(ss, 1)
push!(ss, 2)
push!(ss, 3)

OUT = open("Fx.txt", "w")

for n = 4 : 20000
	m = size(ss, 1)
	tmp = rand(rng, ss)
	(u, v) = edge[tmp]
	push!(edge, (u, n))
	push!(edge, (v, n))
	m += 1
	push!(ss, m)
	m += 1
	push!(ss, m)
	if rlt[n] > 0
		y = getResult(n, edge)
		println(n, " ", y)
		println(OUT, n, " ", y)
	end
end