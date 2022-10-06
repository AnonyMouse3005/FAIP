% Algorithm 6 of Berman (1994)
function [G, E] = alg_berman_mc2(A, s, k, N)
% output: graph G' with (<=k) added edges
arguments
    A (:,:) double  % adj matrix for connected, undirected graph
    s (1,1) double  % s source node
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of agents
end

E = zeros(k, 2);  % list of edges to add
G = graph(logical(A));  % ignore weights of weighted graph
E_G = G.Edges.EndNodes;  % pre-exisiting edges
G = G2SPT(G,s);  % change G to (undirected) SPT
h = ismember(1:length(A),N);  % weight of each node (1 if in N, 0 otherwise)

budget = k;
while budget > 0  % add k edges iteratively (k = Y since uniform weights)
    e = k - budget + 1;
    D = distances(G,s);
    [B, I] = sort(D,'descend');
    testarcs = zeros(1000000,2);
    pos = 1;
    for i = max(B):-1:2  % find all test arcs
        % d(s,j) - d(s,i) >= 0
        arcs = cartprod(I(B<=i), I(B==i));
        idx = ones(1, length(arcs));
        for j = 1:length(arcs)  % remove nonedges (i,j) with i = j and existing edges
            if arcs(j,1) == arcs(j,2) ||...
                    ismember(sort(arcs(j,:)),E_G,'rows') ||...
                    ismember(sort(arcs(j,:)),sort(E,2),'rows') ||...
                    ~ismember(arcs(j,1),N) || ~ismember(arcs(j,2),N)
                idx(j) = 0;
            end
        end
        arcs = arcs(logical(idx),:);
        testarcs(pos:pos+size(arcs,1)-1, :) = arcs;
        pos = pos + size(arcs,1);
    end
    testarcs = testarcs(any(testarcs,2),:);  % remove rows with all zeros
    a = size(testarcs,1);  % # of potential arcs
    
    [u, u_idx] = sort(D(N),'descend');  % M: nodes sorted by nonincreasing order of d(s,j)
    M = N(u_idx);

    if u(1) > u(2)
        improv = alg3_berman(G2SPT(G,s),s,u,M,testarcs);
        [~, bestarc_idx] = max(improv);
        bestarc = testarcs(bestarc_idx, :);  % add links w/ largest savings to G_SPT
        G = addedge(G, bestarc(1), bestarc(2));
        E(e,:) = bestarc;
%         fprintf('iter = %d, add arc (%d,%d)\n', e,bestarc(1),bestarc(2))
        budget = budget - 1;
    else  % u(1) == u(2), step 3
        M_tied = M(u == max(u));  % tied nodes in M for the largest u-value
        t = length(M_tied);
        U = zeros(a, t);
        for ii = 1:a
            aij = testarcs(ii,:);
            for jj = 1:t
                U(ii, jj) = distances(addedge(G,aij(1),aij(2)), s, M_tied(jj));
            end
        end
        zcprime = max(min(U,[],1));
        if max(distances(G,s,N)) == zcprime  % Z_c == Z_c'
            break
        else
            [sol, obj] = IP(G2SPT(G,s), s, testarcs, M_tied, zcprime, h);
            if isempty(sol)
%                 fprintf('cannot solve IP\n')
                break
            end
            if obj > budget  % if R* > Y = remaining budget
                break
            else
                G = addedge(G, sol(:,1), sol(:,2));
                E(e:e+size(sol,1)-1,:) = sol;
                budget = budget - size(sol,1);
                if obj == budget  % if R* = Y
                    break
                end
            end
        end
    end
end

if ~any(E,'all')
    G = graph(logical(A));
else
    E = E(any(E,2),:);
    G = addedge(graph(logical(A)),E(:,1),E(:,2));
end

end

%% Alg 3: return a list of improvement scores for all testarcs
function improv = alg3_berman(SPT, s, u, M, testarcs)

improv = zeros(1,size(testarcs,1));
% step 1
for i = 1:size(testarcs,1)
    arc = testarcs(i,:);
    commonnodes = intersect(shortestpath(SPT,s,arc(1)), shortestpath(SPT,s,arc(2)),'stable');  % nodes common to both paths from s to i and s to j
    [~, idx] = max(distances(SPT,s,commonnodes));
    sprime = commonnodes(idx);
    path = shortestpath(SPT,sprime,arc(2));  % path p' from s' to j
    S = zeros(size(path));
    % step 2
    for j = 1:length(path)  % j equiv to k in paper
        Dk = distances(SPT,s,path(j)) - (distances(SPT,s,arc(1)) + 1 + distances(SPT,arc(2),path(j)));
        if Dk > 0
            S(j) = Dk;
        end
    end
    ubar = zeros(1, length(M));  % stores ubar values
    path_noXprime = path;
    if sprime == s  % ignore s in p' if sprime = s
        path_noXprime = setdiff(path,s,'stable');
    end
    % step 3
    for alpha = 1:length(M)
        if ~any(ismember(shortestpath(SPT,s,M(alpha)), path_noXprime))  % path from s to m_alpha does not include nodes in p'
            ubar(alpha) = u(alpha);
            break
        else
            for k = length(path):-1:1  % find kprime
                if ismember(path(k), shortestpath(SPT,s,M(alpha)))
                    jprime = k;
                end
            end
            if S(jprime) == 0
                ubar(alpha) = u(alpha);
                break
            else
                % step 4
                ubar(alpha) = min([u(alpha), distances(SPT,s,M(alpha))-S(jprime)]);
                if ubar(alpha) == u(alpha)
                    break
                end
            end
        end
    end
    % step 5
    improv(i) = u(1) - max(ubar);
end

end

%% IP
function [E, obj] = IP(SPT, s, testarcs, M, zcprime, h)

m = length(M);
a = size(testarcs,1);
A = zeros(m, a);  % LHS of inequality constraints (in form a*x <= b)
for i = 1:a
    aij = testarcs(i,:);
    for j = 1:m
        mk = M(j);
        if h(mk)*distances(addedge(SPT,aij(1),aij(2)), s, mk) < zcprime
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
E = testarcs(x==1,:);

end

%% Convert general graph to undirected SPT
    function SPT = G2SPT(G, s)
        SPT = shortestpathtree(G,s);
        SPT = graph(SPT.Edges);
    end
