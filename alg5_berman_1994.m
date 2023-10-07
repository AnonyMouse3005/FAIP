% Algorithm 5 of Berman (1994)
function [G, E] = alg5_berman_1994(A, s, k, N, zcprime)
% output: graph G' with (<=k) added edges
arguments
    A (:,:) double  % adj matrix for connected, undirected graph
    s (1,1) double  % pre-located facility
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of agents
    zcprime (1,1) double = 5  % initial target MAC
end

G = graph(logical(A));  % ignore weights of weighted graph
V = 1:length(A);
nonEdges = setdiff(nchoosek(N,2),G.Edges.EndNodes,'rows');  % set of non-edges in G (a x 2)
a = length(nonEdges);  % number of nonedges

G_SPT = shortestpathtree(G,s);  % change G to (undirected) SPT
G_SPT = graph(G_SPT.Edges);
d = distances(G_SPT,s);
h = ismember(1:length(A),N);  % weight of each node (1 if in N, 0 otherwise)

count = 0;
obj = [];
zcp_history = [zcprime, zeros(1,1000)];
sol = 0;
while isempty(obj) || all(obj ~= k)
    if ~issorted(zcp_history(any(zcp_history,1)),'monotonic')  % break if Z'c go back and forth
%         fprintf('Cannot find solution!\n')
        break
    end
    if count ~= 0
        if isempty(obj) || obj > k
            zcprime = zcprime + 1;
            zcp_history(count+1) = zcprime;
        elseif ~isempty(obj) && obj < k
            zcprime = zcprime - 1;
            zcp_history(count+1) = zcprime;
        end
    end
    M = V(d>=zcprime);  % critical nodes (distances to s >= target max cost)
    m = length(M);  % number of critical nodes 
    A = zeros(m, a);  % LHS of inequality constraints (in form a*x <= b)
    for i = 1:a
        e = nonEdges(i,:);
        for j = 1:m
            mk = M(j);
            if h(mk)*distances(addedge(G_SPT,e(1),e(2)), s, mk) < zcprime
                A(j,i) = -1;
            end
        end
    end
    f = ones(a,1);  % objective function
    intcon = 1:a;  % integer conditions
    b = -1*ones(m,1);  % RHS of inequality constraints
    lb = zeros(a,1);  % lower bound and upper bound enforce all variables are binary
    ub = ones(a,1);
    options = optimoptions('intlinprog','Display','none');
    [x,obj,~,~] = intlinprog(f,intcon,A,b,[],[],lb,ub,[],options);
    if any(x)
        sol = x;
    end
    count = count+1;
end

E = nonEdges(sol==1,:);
G = addedge(G,E(:,1),E(:,2));

end
