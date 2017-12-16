function [graph_exists, A, D] = checkout_adjmx(graph_name, data_dir)
% function [graph_exists, A, D] = checkout_adjmx(graph_name)
%
% Using fixed data directory structure, get an adjacency matrix for a
% particular graph.
%
% If adjacency matrix already exists in the adjacency_matrices directory,
% it will be loaded and returned as A (along with degree matrix D). If not,
% A and D are set as scalar zero, and graph_exists is returned false.

adj_mx_dir = '/adjacency_matrices';
thepath = strcat(data_dir, adj_mx_dir, '/', graph_name, '.mat');

% adjacency_mxs_struct = dir(strcat(data_dir, '/adjacency_matrices'));
% adjacency_mxs = {adjacency_mxs_struct.name};
% adjacency_mx_fns = {};
% for i=1:length(adjacency_mxs)
%   if ~isempty(regexp(adjacency_mxs{i}, '.*.graph'))
%     [~, adjacency_mx_fns{end+1}, ~] = fileparts(adjacency_mxs{i});
%   end
% end
% 
% graph_matches = strcmp(graph_name, adjacency_mx_fns);
% graph_exists = any(graph_matches);

graph_exists = exist(thepath, 'file');

if graph_exists
  % Load it
  data = load(thepath);
  A = data.A;
  D = data.D;
else
  A = 0;
  D = 0;
end