include("Graph.jl")
include("GraphMatrix.jl")
include("Algorithm.jl")
include("SetCreation.jl")
include("Experiment.jl")

# Read Graph
buf = split(ARGS[1], ',')
fileName = string(buf[1], ".txt")
networkType = "unweighted"
if size(buf, 1) > 1
	networkType = buf[2]
end
# Find LLC
G0 = readGraph(fileName, networkType)
G = getLLC(G0)
# Read n2,k,Methold
buf2 = split(ARGS[2], ',')
n2 = parse(Int32, buf2[1])
k = parse(Int32, buf2[2])
Methold = "r"
if size(buf2, 1) > 2
	Methold = buf2[3]
end
# Create Q, EQ
Q = createQ(G, n2, Methold)
EQ = createEQ(G, Q)
# Do the Experiment
algs = split(ARGS[3], ',')
expNum = size(algs, 1)
OUT = open("Result.txt", "a")
println(OUT, buf[1], " ", G.n, " ", G.m)
println(OUT, n2, " ", k)
println(OUT, expNum)
for i = 1 : expNum
	S, T = runExperiment(G, Q, EQ, k, algs[i])
	if (expNum > 1)
		println(OUT, algs[i], " ", T, " ", EffectiveResistance(G, Q, S))
	else
		println(OUT, algs[i], " ", T)
	end
end
