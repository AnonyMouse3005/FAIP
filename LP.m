% Implementation of Charikar's LP-based algorithm for k-median w/ penalties
function G = LP(A, s, k, N)
% output: resultant graph G after adding k edges
arguments
    A (:,:) double  % adj matrix (sparse) for connected, undirected graph
    s (1,1) double  % source node
    k (1,1) double  % k facilities to open/k edges to be added
    N (1,:) double  % node subgroup i.e., clients
end

G = graph(logical(A));
if k >= length(N)
    return;
end

V = 1:length(A);  % facilities F=V
v = length(V);
n = length(N);
c = distances(G,V,N);  % cost matrix
p = zeros(1,n);  % penalty array
for j = 1:n
    p(j) = distances(G,N(j),s);
    for i = 1:v
        if i ~= s && N(j) ~= s
            c(i,j) = c(i,j) + 1;  % add 1 to cost if neither i nor j is s
        end
    end
end

obj = [reshape(c',[],1); p'; zeros(v,1)];  % obj function
intcon = 1:length(obj);  % integer conditions
ub = ones(size(obj));  % upper bounds (all 1)
lb = zeros(size(obj));  % lower bounds (all 0)
A = zeros(v*n+n,length(obj));  % constraint matrix
b = zeros(v*n+n,1);
for a = 1:v  % first inequality constraint
    for aa = 1:n
        A((a-1)*n+aa,(a-1)*n+aa) = 1; A((a-1)*n+aa, v*n+n+a) = -1;
    end
end
for a = 1:n  % second inequality constraint
    A(v*n+a,v*n+a) = -1;
    b(v*n+a) = -1;
    for aa = a:n:v*n
        A(v*n+a,aa) = -1;
    end
end
Aeq = [zeros(1,v*n+n),ones(1,v)];  % linear equality constraint
beq = k;
options = optimoptions('intlinprog','Display','none');
x = intlinprog(obj,intcon,A,b,Aeq,beq,lb,ub,options);

% xij = reshape(x(1:v*n),[n,v])';
% r = x(v*n+1:v*n+n)';
y = x(v*n+n+1:end)';

u = V(y>0);  % solution, k facilties opened

G = addedge(G, ones(1,length(u))*s, u);

end