function [A, D] = graph_to_adjmx(graph_name, data_dir)
% function [A, D] = graph_to_adjmx(graph_name, data_dir)
%
% Given graph name and data directory, load the graph from CHACO/METIS format
% and create the adjacency matrix using matlab's sparse matrix format.

graph_dir = '/graphs';
thepath = strcat(data_dir, graph_dir, '/', graph_name, '.graph');

f = fopen(thepath, 'r+', 'native', 'US-ASCII');

tstart = tic;

if f > 0
  r = fseek(f, 0, 'bof');
  if r < 0
    fclose(f);
    error('  Seeking in file failed.');
  end
  info = split(fgetl(f));
  if length(info) == 2
    n = str2double(info(1));
    m = str2double(info(2));
    
    start_idx = 1;
    stride = 1;
  elseif length(info) == 3
    warning('  Graph has weights; may not work correctly');
    % Code below just ignores the weights anyways.
    n = str2double(info(1));
    m = str2double(info(2));
    if length(info(3)) == 1
      weightc1 = 0;
      weightc2 = str2double(info{3}(1));
    else
      weightc1 = str2double(info{3}(1));
      weightc2 = str2double(info{3}(2));
    end
    
    if weightc1 == 0
      start_idx = 1;
    else
      start_idx = 2;
    end
    if weightc2 == 0
      stride = 1;
    else
      stride = 2;
    end
  else
    fclose(f);
    error('  Info string (first line of file) did not have length 2 or 3; unrecognized format');
  end
  
  % Create sparse matrix index data
  row_data = zeros(2*m, 1);
  col_data = zeros(2*m, 1);
  
  
  % Read graph data
  vertex_id = 1;
  index_ptr = 1;
  while 1
    line = fgetl(f);
    line = line(2:end); % get rid of whitespace offsetting
    if ~ischar(line)
      % break if eof -- graph has been read
      break;
    end
    data = split(line);
    neighbor_idxs = start_idx:stride:length(data);
    neighbors = str2double(data(neighbor_idxs));
    row_data(index_ptr:index_ptr+length(neighbors)-1) = neighbors;
    col_data(index_ptr:index_ptr+length(neighbors)-1) = vertex_id*ones(length(neighbors),1);
    
    % status
    if mod(vertex_id, 5e3) == 0
      fprintf('  (time %f) processed vertex %d of %d\n', toc(tstart), vertex_id, n);
    end
    % updates
    vertex_id = vertex_id + 1;
    index_ptr = index_ptr + length(neighbors);
  end
%   keyboard
  fclose(f);
else
  fprintf('  Failed to open file %s with error %d\n', thepath, f);
end

fprintf('  (time %f) finished reading graph data\n', toc(tstart));

A = sparse(row_data, col_data, ones(2*m, 1), n, n);
D = spdiags(sum(A, 2), 0, n, n);