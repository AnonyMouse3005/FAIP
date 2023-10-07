% Implementation of FFT and FFT_V
function G = FFT(A, s, k, N)
% output: resultant network G after adding k edges
arguments
    A (:,:) double  % adj matrix (sparse) for connected, undirected network
    s (1,1) double  % pre-located facility
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of agents; when N = V, we have FFT_V
end

u = zeros(1, k);
G = graph(logical(A));
N = sort(N);
D = distances(G,s,N);  % distance from s to all nodes in N

i = 1;
while i <= k
    [~, idx] = max(D);
    u(i) = N(idx);
    d = distances(G,u(i),N);
    for j = 1:length(N)
        D(j) = min(D(j), d(j));
    end
    i = i+1;
end

G = addedge(G, ones(1,length(u))*s, u);

end