% Implementation of Charikar's LP-based algorithm for k-median w/ penalties
function G = LP(A, s, k, N)
% output: resultant network G after adding k edges
arguments
    A (:,:) double  % adj matrix (sparse) for connected, undirected network
    s (1,1) double  % pre-located facility
    k (1,1) double  % k extra facilities to open/k edges to be added
    N (1,:) double  % set of agents
end

G = graph(logical(A));
if k >= length(N)
    G = addedge(G, ones(1,length(N))*s, N);
    return;
end

v = length(A);
V = 1:v;  % facilities F=V
n = length(N);
c = distances(G,V,N);  % cost matrix (v x n)
p = zeros(1,n);  % penalty array
for j = 1:n  % define cost matrix and penalty array
    p(j) = distances(G,N(j),s);
    for i = 1:v
        if i ~= s && N(j) ~= s
            c(i,j) = c(i,j) + 1;  % add 1 to cost if neither i nor j is s
        end
    end
end

z_min = 0;
z_max = n;
z1 = 999; k1 = 0; k2 = n;  % '2' for current results, '1' for previous results
z2 = z_min + floor((z_max-z_min)/2);
counter = 0;  % # of facility location w/ penalty algorithm calls
while true
    if k2 ~= k && abs(z1-z2) > 1  % perform binary search on z2
        [I2, flp2, k2] = FLP_JV(v,n,z2,c,p,N);
        counter = counter+1;
        if k2 < k
            z_max = z2 - 1;
        else  % k2 > k
            z_min = z2;
        end
        if z_max - z_min > 2
            [z1, I1, flp1, k1] = deal(z2, I2, flp2, k2);
            z2 = z_min + floor((z_max-z_min)/2);
        else  % solution is close
            if z_min == z_max
                range = z_max;
            else
                range = z_min+1:z_max;
            end
            for z2 = range
                if abs(z1-z2) <= 1 && ((k2 >= k && k1 <= k) || (k2 <= k && k1 >= k))
                    break
                end
                [I1, flp1, k1] = deal(I2, flp2, k2);
                [I2, flp2, k2] = FLP_JV(v,n,z2,c,p,N);
                z1 = z2;
                counter = counter+1;
            end
        end
    elseif k2 ~= k && abs(z1-z2) <= 1  % do augmentation step
        if length(I1) < length(I2)
            I_final = Aug(I1,flp1,I2,flp2,k);
        else
            I_final = Aug(I2,flp2,I1,flp1,k);
        end
        break
    else  % k2 == k
        I_final = I2;
        break
    end
end
% fprintf('Charikar Algorithm for n=%d, k=%d done after %d iterations!\n',n,k,counter)
% disp(I_final)
G = addedge(G, ones(1,length(I_final))*s, I_final);

end

%% Algorithm for facility location with penalties by Jain-Vazirani (1999)
function [I_open, flp, k] = FLP_JV(v,n,z,c,p,N)
%     fprintf('Start FLP instance for z=%d\n',z)
    flp = FLP(v,n,z,c,p);
    
    while ~all(flp.connected)  % terminate when all clients are connected/assigned conwits
        J = find(~flp.connected);
        if isequal(flp.alpha(J), flp.penalty(J))  % stop early if alpha cannot further increase (bounded by p)
%             fprintf('Phase 1 stopped early: %d unconnected clients\n',n-sum(flp.connected))
            break
        end
        flp.growAlpha();
        flp.t = flp.t + 1;
    end
%     fprintf('Phase 1 done! z=%d: %d paid facilities\n',z,sum(flp.temp_open))
    A_tight = zeros(v,v);
    A_tight(:,N) = logical(flp.beta);
    G_tight = simplify(graph(logical(A_tight + A_tight')));  % create subgraph from saturated edges
    I_paid = find(isfinite(flp.t_paid));
    [~, t_idx] = sort(flp.t_paid(I_paid));
    I_paid = I_paid(t_idx);  % sort facilities by paid-for time
    I_open = zeros(size(I_paid));
    count = 1;
    while ~isempty(I_paid)
        I_open(count) = I_paid(1);
        I_near = nearest(G_tight,I_paid(1),2);
        I_paid = setdiff(I_paid, I_near);  % delete facilities reachable from i by <= 2 tight edges
        I_paid(1) = [];
        count = count + 1;
    end
    I_open(I_open==0) = [];
%     kprime = length(I_open);
%     fprintf('Phase 2 done! z=%d: %d open facilities\n',z,kprime)
    k = length(I_open);
end

%% Augmentation step
function I_final = Aug(I1,flp1,I2,flp2,k)  % I1 small sol.; I2 large sol.
    I_extra = setdiff(I2,I1);
    I_final = zeros(1,k);
    I_final(1:length(I1)) = I1;  % adding centers from large sol to small sol
    i = length(I1)+1;
    add = true;
    while ~isempty(I_extra) && nnz(I_final) ~= k
        stop = false;
        for iprime = 1:length(I1)
            if stop
                break
            end
            for j = 1:length(flp1.alpha)
                if flp1.beta(I_extra(1),j)>0 && flp2.beta(I_extra(1),j)>0 ...
                        && flp1.beta(I1(iprime),j)>0 && flp2.beta(I1(iprime),j)>0
                    add = false;
                    stop = true;
                    break
                end
            end
        end
        if add
            I_final(i) = I_extra(1);
        end
        I_extra(1) = [];
        i = i + 1;
    end
    I_final(I_final==0) = [];  % in case isempty(I_extra)
end