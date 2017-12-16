function [S] = compute_sparsest_cut_tree(A, D)
% function [S] = compute_sparsest_cut_tree(A, D)
%
% Given adjacency and degree matrix for a tree A,D, find sparsest cut S
%
% If A/D are not from a tree this will not work correctly.

