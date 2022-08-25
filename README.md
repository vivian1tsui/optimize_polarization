# Minimizing Polarization in Noisy Leader-Follower Opinion Dynamics

This code accompanies the paper "Minimizing Polarization in Noisy Leader-Follower Opinion Dynamics" submitted to IEEE Transactions on Knowledge and Data Engineering. The code is written in *Julia v1.3.1*.

## Datasets

We select a large set of real-world social networks from [Koblenz Network Collection](http://konect.uni-koblenz.de/) and [SNAP](https://snap.stanford.edu). All networks tested are regarded as undirected networks. For those networks that are disconnected originally, we perform our experiments on their largest connected components.

## Usage

### Network Input

The input of a network consists of the edges in the network. Each line of the input file represents a undirected unweighted edge in the network, which is specified as the format "source_node target_node". Some example files have been put in the `Data/` repository.

All input files should be put in the `Data/` and have a filename ending with '.txt'. The algorithms will run on these files in lexicographical order of the file names.

### Run

Execute command `julia Test.jl ARGS[1] ARGS[2] ARGS[3]`, where
the directory where the edges list files are.

- `ARGS[1]` : networkName, networkType
  `ARGS[2]` : n2, k, selectMethod
  `ARGS[3]` : Algorithm[1], Algorithm[2], Algorithm[3],...
- `networkType` = {"weighted", "unweighted"}
  `n2` = |Q|
  `selectMethold` = {'r', 'l', 's'}  #Choose Q Method == 'r'(random) | 'l'(largest degree) | 's'(smallest degree)
  `Algorithm[i]` = {"Optimal", "Random", "topCent", "topDegree", "Exact", "Approx:eps"}
- The result will be printed to both console and file `result.txt`.

### Examples

```bash
# Run the  approximate algorithms with $\epsilon$=0.5, the graph has |Q|=10 leaders and randomly choosed k=20 edges adding to Q, 
$ julia Test.jl Data/Adjnoun 10,20 Approx:0.5

```

