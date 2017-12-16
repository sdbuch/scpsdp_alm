function [D] = pdistmx(X)
%% [D] = pdistmx(X)
% Return pairwise distances given a data matrix X
% Inputs:
% - X: d by n matrix of d-dimensional data points
% Outputs:
% - D: n by n matrix of pairwise distances, where:
%      d_ij = norm(X(:,i) - X(:,j), 2) for i,j = 1, ..., n

n = size(X,2);
d = size(X,1);
for i = 1:n
  for j = 1:n
    D(i,j) = norm(X(:,i) - X(:,j), 2);%.^2 * 1/2;
  end
end

return