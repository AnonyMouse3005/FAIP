% Implementation of Local Search Algorithm for k-median w/ penalties
function G = LS(A, s, k, N, q, delta)
% output: resultant network G after adding k edges
arguments
    A (:,:) double  % adj matrix (sparse) for unweighted, undirected network
    s (1,1) double  % source node
    k (1,1) double  % k number of edges to be added
    N (1,:) double  % node subgroup
    q (1,1) double = 1  % number of edges to be swapped in each iteration
    delta (1,1) double = 0.95  % between 0 and 1, the closer to 1 the better the approximation
end

seed = 1;  % seed for generating Eprime arbitrarily
v = size(A,1);
G = graph(logical(A));
s_nebors = neighbors(G,s)';  % neighbors of s
s_nnebors = setdiff(1:v,s_nebors);  % non-neighbors of s
uCand = s_nnebors;
rng(seed)
u = uCand(randsample(length(uCand),k));  % arbitrary solution set of size k (endnodes of edges in E')
m = length(uCand) - k;

improvable = 1;
u1_comb = nchoosek(1:k,q);  % set of all combinations of indices of (E')'s edges
u2_comb = nchoosek(1:m,q);  % set of all combinations of indices of (\bar{E}^V_S - E')'s edges
while 1
    udiff = setdiff(uCand, u);  % endnodes of edges in \bar{E}^V_S - E'
    better = 0;  % flag to break both for loops i.e., found better solution set
    for i = 1:size(u1_comb,1)
        u1 = u(u1_comb(i));  % subset of u of size q
        for j = 1:size(u2_comb,1)
            u2 = udiff(u2_comb(j));  % subset of udiff of size q
            utemp = [setdiff(u, u1), u2];  % [(E'-E1) \cup E2]
            TNSEtemp = calcTNSE(addedge(G, ones(1,length(utemp))*s, utemp), s, N);
            TNSEprime = calcTNSE(addedge(G, ones(1,length(u))*s, u), s, N);
            if TNSEtemp <= delta*TNSEprime  % larger delta (up to 1), slower run, better approx
                better = 1;
                break;
            end
            if i == size(u1_comb,1) && j == size(u2_comb,1)
                improvable = 0;  % solution cannot be improved further
            end
        end
        if better == 1
            break;
        end
    end
    if improvable == 0
        break;
    end
    u = utemp;
end

G = addedge(G, ones(1,length(u))*s, u);

end
