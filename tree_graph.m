% make a tree graph's adjacency matrix
% Just use a minimum spanning tree from a random graph
n = 10;
p = 0.5; %erdos-renyi parameter

edges = rand(n,n);
edges = tril(edges) - diag(diag(edges));
edge_mask = edges > 1-p;
edge_mask = edge_mask + edge_mask';
A = double(edge_mask);

G = graph(A);
T = minspantree(G);

AA = T.adjacency;

imagesc(A);
