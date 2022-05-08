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

if ARGS[2] == "createQ"
	# Create the Set Q,EQ
	buf2 = split(ARGS[3], ',')
	n2 = parse(Int32, buf2[1])
	Methold = "r"
	if size(buf2, 1) > 1
		Methold = buf2[2]
	end
	Q = createQ(G, n2, Methold)
else
	buf2 = split(ARGS[2], ',')
	k = parse(Int32, buf2[1])
	algs = split(ARGS[3], ',')
	Q = readQ(G)
	n2 = size(Q, 1)
	EQ = createEQ(G, Q)
	# Do the Experiment
	expNum = size(algs, 1)
	ALLT = Array{Float64}(expNum)
	ALLS = Array{Array{Tuple{Int32, Int32, Float64}, 1}}(expNum)
	for i = 1 : expNum
		ALLS[i], ALLT[i] = runExperiment(G, Q, EQ, k, algs[i])
	end
	# Output the result
	OUT = open("Result.txt", "a")
	println(OUT, buf[1], " ", G.n, " ", G.m)
	println(OUT, n2, " ", k)
	println(OUT, expNum)
	printResult(G, Q, OUT, algs, ALLS, ALLT)
	println(OUT)
end
