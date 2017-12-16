function [success] = checkin_adjmx(graph_name, data_dir, A, D)
% function [success] = checkin_adjmx(graph_name, data_dir, A, D)
%
% Using fixed data directory structure, save an adjacency matrix for a
% particular graph. (created with graph_to_adjmx.m)
%
% If adjacency matrix already exists in the adjacency_matrices directory,
% success is set to -1. Other failures set success to 0. Else the file is
% written (.mat file) and success is set to 1);

adj_mx_dir = '/adjacency_matrices';
thepath = strcat(data_dir, adj_mx_dir, '/', graph_name, '.mat');

graph_exists = exist(thepath, 'file');
if graph_exists
  success = -1;
  return;
else
  save(thepath, 'A', 'D');
  success = 1;
  return;
end