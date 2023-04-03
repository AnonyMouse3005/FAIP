% Implementation of R-Prop and R-Uniform
function u = RandF(A, s, k, N, w, replace, nrep)
% output: matrix u containing candidate nodes to be connected to s, sorted by order of adding,
% each row is one run
arguments
    A (:,:) double  % adj matrix (sparse) for connected, undirected network
    s (1,1) double  % pre-located facility
    k (1,1) double  % k edges to be added
    N (1,:) double  % set of clients
    w (1,1) double = 1  % sampling with weighted probs proportional to agents' relative distances from s (1) or uniformly (0)
    replace (1,1) double = 0  % sampling with (1) or without (0) replacement
    nrep (1,1) double = 10  % number of runs
end

if k > length(N)
    disp('k capped at |N|. Will use k = |N| instead.')
    k = length(N);
end

u = zeros(nrep, k);
G = graph(A);
N = sort(N);
D = distances(G,s,N);  % distance from s to all nodes in N

if replace  % sampling w/ replacement
    if w  % sample with weighted probabilities
        for r = 1:nrep
            rng(r)
            u(r,:) = randsample(N,k,true,D/sum(D));
        end
    else  % sample uniformly
        for r = 1:nrep
            rng(r)
            u(r,:) = randsample(N,k,true);
        end
    end
else  % sampling w/o replacement
    if w  % sample with weighted probabilities (NOT USED IN PAPER)
        for r = 1:nrep
            % Note: randsampleFS also selects those with 0 weights (after
            % all agents w/ positive weights have been selected)
            % Therefore, preprocess x to contain only agents w/ positive weights
            wgts = D/sum(D);
            p = N(wgts>0);
            if k >= length(p)
                u(r,:) = p;
                continue
            end
            rng(r)
            y = randsampleFS(length(p),k,wgts(wgts>0));
            u(r,:) = p(y);
        end
    else  % sample uniformly
        for r = 1:nrep
            rng(r)
            u(r,:) = randsample(N,k);
        end
    end
end

end