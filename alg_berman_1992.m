% First heuristic from page 12 of Berman (1992)
function [G, E] = alg_berman_sc(A, s, k, N)
% output: graph G' with k added edges
arguments
    A (:,:) double  % adj matrix for connected, undirected graph
    s (1,1) double  % s source node
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of agents
end

E = zeros(k, 2);  % list of edges to add
G = graph(logical(A));  % ignore weights of weighted graph

G = G2SPT(G,s);  % change G to (undirected) SPT
for j = 1:k  % add k edges iteratively (k = Y since uniform weights)
    D = distances(G,s);
    [B, I] = sort(D,'descend');
    testarcs = zeros(1000000,2);
    pos = 1;
    for i = max(B):-1:2  % find all test arcs
        % 1st test: d(s,j) - d(s,i) > l(i,j) = 1 (unweighted graph)
        arcs = cartprod(I(B<i-1), I(B==i));
        testarcs(pos:pos+size(arcs,1)-1, :) = arcs;
        pos = pos + size(arcs,1);
    end
    % 2nd test is not needed for unweighted graph
    testarcs = testarcs(any(testarcs,2),:);  % remove rows with all zeros
    
    improv = cutset(G2SPT(G,s),s,testarcs,N);
    [~, bestarc_idx] = max(improv);
    bestarc = testarcs(bestarc_idx, :);
    G = addedge(G, bestarc(1), bestarc(2));

    E(j,:) = bestarc;
end

G = addedge(G,E(:,1),E(:,2));

end

%% cutset method
function improv = cutset(SPT, s, testarcs, N)  % return improvement scores for all testarcs

improv = zeros(1, length(testarcs));
for i = 1:length(testarcs)  % for each test arc
%     fprintf("i=%d\n",i)
    arc = testarcs(i,:);
    [path, d] = shortestpath(SPT,s,arc(2));
    for k = 1:d
        if k == d
            improv(i) = 0;
            break
        end
        sav = distances(SPT,s,path(k+1)) - (distances(SPT,s,arc(1)) + 1 + distances(SPT,arc(2),path(k+1)));
        if sav > 0
%             fprintf('[%s]\n', join(string(path), ','));
%             fprintf("arc (%d,%d), sav=%d = %d-[%d+1+%d], W_k=%d\n",arc(1),arc(2),sav, ...
%                 distances(SPT,s,path(k+1)), distances(SPT,s,arc(1)), distances(SPT,arc(2),path(k+1)),path(k+1))
            bins = conncomp(SPT.rmedge(path(k),path(k+1)));
            nodes = 1:SPT.numnodes;
            nodes = nodes(bins==bins(path(k+1)));  % nodes in subtree that contains node W_k
%             fprintf('[%s]\n', join(string(nodes), ','));
            for v = nodes
                if ismember(v, N)
                    improv(i) = improv(i) + distances(SPT,s,v) - (distances(SPT,s,arc(1)) + 1 + distances(SPT,arc(2),v));
                end
            end
%             fprintf('improv=%d\n',improv(i))
            break
        end
    end
end

end

%% Convert general graph to undirected SPT
    function SPT = G2SPT(G, s)
        SPT = shortestpathtree(G,s);
        SPT = graph(SPT.Edges);
    end
