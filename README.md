# FAIP

This repo contains MATLAB source codes implemented for all proposed algorithms and benchmarks for facility accessibility improvement problems (FAIPs), along with datasets and all resulting figures from the experiments.


### Algorithms

(all take sparse adjacency matrix `A` of network $G=(V,E)$, source node `s`, number of edges to add `k`, and set of clients `N`$\subseteq V$ as inputs)
- `FFT.m`: Farthest-First Traversal.
- `k_Im.m`: Highest $k$ Importance. Additionally requires specifying the name of a centrality measure ("random", "highdegree", "lowdegree", "betweenness", "pagerank", "eigenvector", "clusteringcoefficient", "closeness", "eccentricity").
- `RandF.m`: Randomized Framework. Optionally specify `w`, `replace` , `nrep` parameters to run sampling with weighted probabilities or uniformly, with or without replacement, and with how many number of runs (more details in documentation).
- `LS.m`: Local Search Algorithm in [5]. `q` and `delta` parameters to specify the number of edges to be swapped in each iteration and the quality of the solution.
- `LP.m`: Linear Programming-based Algorithm using primal-dual schema in [3]. The object for storing the facility location with penalty instance is called from `FLP.m`.
- `alg_berman_1992.m`: First heuristics in page 12 of [2].
- `alg5_berman_1994.m`: Algorithm 5 in [1]. `zcprime` parameter to intialize a target MNSE of choice.
- `alg6_berman_1994.m`: Algorithm 6 in [1].

### Helper functions
- `calcMAC.m`: Compute maximum accessibility cost of a network `G` (graph object) with source `s` and set of clients `N`.
- `calcTAC.m`: Compute total accessibility cost of a network `G` (graph object) with source `s` and set of clients `N`.
- `lib/ClusteringCoefficient.m`: Compute (local) clustering coefficient of network $G$, adapted from [URL](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45734/versions/1/previews/cnm/avgClusteringCoefficient.m/index.html).

### Datasets
- `datasets/dt_oregon.mat`: AS peering networks, available at [4].
- `datasets/p2p_connected.mat`: Internet peer-to-peer networks, available at [4].
- `datasets/social_small.mat`: Small social networks, available at [6].
- `datasets/erdos_renyi_200.mat`: Small Erdos-Renyi random networks (200 nodes).
- `datasets/erdos_renyi_3000.mat`: Large Erdos-Renyi random networks (3000 nodes).

<br/><br/>


### References

[1] O. Berman, D. I. Ingco, and A. Odoni, _Improving the location of minimax facilities through network modification_, Networks, 24 (1994), pp. 31–41.

[2] O. Berman, D. I. Ingco, and A. R. Odoni, _Improving the location of minisum facilities through network modification_, Ann. Oper. Res., 40 (1992), pp. 1–16.

[3] M. Charikar, S. Khuller, D. M. Mount, and G. Narasimhan, _Algorithms for facility location problems with outliers_, in SODA, vol. 1, 2001, pp. 642–651.

[4] J. Leskovec and A. Krevl, _SNAP Datasets: Stanford large network dataset collection_. [URL](http://snap.stanford.edu/data/), June 2014.

[5] A. Meyerson and B. Tagiku, _Minimizing average shortest path distances via shortcut edge addition_, in Approx/Random, Springer, 2009, pp. 272–285.

[6] R. A. Rossi and N. K. Ahmed, _The network data repository with interactive graph analytics and visualization_, in AAAI, 2015. [URL](https://networkrepository.com/).
