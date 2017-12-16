function [v] = cell_helper(cell_arr, gram_mx)
% function [v] = cell_helper(cell_arr, gram_mx)
%
% Helper function for use with BM solver
%
% cellfun(@(x3) ( sum( vec( x3 .* (Y*Y')) ) ), triangle_cstr) does not work
% when triangle_cstr contains sparse matrices. So just do the same thing
% here with a function.

v = zeros(size(cell_arr, 1), 1);
for k = 1:length(v)
  v(k) = sum(vec(cell_arr{k} .* gram_mx));
end