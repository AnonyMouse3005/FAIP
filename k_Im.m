% Implementation of k-Im and k-Im_V
function G = k_Im(A, s, k, N, m)
% output: resultant network G after adding k edges
arguments
    A (:,:) double  % adj matrix (sparse) for connected, undirected network
    s (1,1) double  % pre-located facility
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of clients; when N = V, we have k-Im_V
    m (1,1) string  % a string in ["random", "highdegree", "lowdegree",...
%         "betweenness", "pagerank", "eigenvector",...
%         "clusteringcoefficient", "closeness", "eccentricity"]
end

G = graph(logical(A));
v = length(A);

I = zeros(v,1);  % node importance (column) array
switch m
    case "clusteringcoefficient"
        I = ClusteringCoefficient(A);
    case "random"
        nrep = 10;
        for j = 1:nrep
            rng(j)
            I = I + randperm(v).';
        end
        I = I / nrep;  % avg over nrep runs
    case "highdegree"
        I = centrality(G, "degree");
    case "lowdegree"
        I = 1./centrality(G, "degree");
    case "eccentricity"
        I = 1./max(distances(G),[],2);
    otherwise
        I = centrality(G, m);
end
I = I.';

u = I(sort(N));  % only keep importance measures from candidate nodes
[~, sortIdx] = sort(u,'descend');
u = N(sortIdx);  % rank the candidate nodes
u = u(1:k);

G = addedge(G, ones(1,length(u))*s, u);

end