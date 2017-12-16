%% script for running sparsest cut on a graph
clear variables; close all;
%% params
% directory structure
% (note some directory structure is set in checkout/checkin_adjmx.m files)
root_dir = '~/projects/matlab/bm_semimetric';
scripts_dir = root_dir;
data_dir = strcat(root_dir, '/data');

graph_name = 'tree40';
semimetric_mode = 'sdp';

tstart = tic;
%% Graph loading
% Check if we already have an adjacency matrix for this graph
[graph_exists, A, D] = checkout_adjmx(graph_name, data_dir);

if graph_exists
  % Adjacency matrix exists, and we loaded it
  fprintf('(time %f) Adjacency matrix for graph %s loaded.\n', toc(tstart), graph_name);
else
  % No adjacency matrix exists. Create it from the .graph file
  fprintf('(time %f) Adjacency matrix for graph %s doesn''t exist; creating it.\n',...
    toc(tstart), graph_name);
  [A, D] = graph_to_adjmx(graph_name, data_dir);
  checkin_adjmx(graph_name, data_dir, A, D);
end

%% Semimetric learning phase
switch semimetric_mode
  case 'lp'
    % LINEAR PROGRAM SOLVER
    % Use MATLAB's linprog solver.
    fprintf('(time %f) Using the LP solver.\n', toc(tstart));
    fprintf('(time %f) Creating constraint matrices.\n', toc(tstart));
    n = size(A,1);
    % there are nchoosek(n, 2) decision variables in the problem
    triangle_idxs = zeros(3*3, nchoosek(n,3));
    coord_map = @(i, j)((i-1)*(2*n - i)/2 + (j - i)); % return an idx into [nchoosek(n,2)] given a pair (i, j), i<j
    counter = 1;
    %     keyboard
    for i = 1:n
      for j = i+1:n
        for k = j+1:n
          % i < j < k
          % so calls to coord_map need to respect this ordering
          %           v1 = sparse(...
          %             [coord_map(i,j); coord_map(i, k); coord_map(j, k)],...
          %             ones(3,1), [1; 1; -1], nchoosek(n, 2), 1);
          %           v2 = sparse(...
          %             [coord_map(i,k); coord_map(j, k); coord_map(i, j)],...
          %             ones(3,1), [1; 1; -1], nchoosek(n, 2), 1);
          %           v3 = sparse(...
          %             [coord_map(i,j); coord_map(j, k); coord_map(i, k)],...
          %             ones(3,1), [1; 1; -1], nchoosek(n, 2), 1);
          v1 = [coord_map(i,j); coord_map(i, k); coord_map(j, k)];
          v2 = [coord_map(i,k); coord_map(j, k); coord_map(i, j)];
          v3 = [coord_map(i,j); coord_map(j, k); coord_map(i, k)];
          %           keyboard
          triangle_idxs(1:3, counter) = v1;
          triangle_idxs(4:6, counter) = v2;
          triangle_idxs(7:9, counter) = v3;
          counter = counter + 1;
          %           if mod(counter, floor(nchoosek(n,3)/10)) == 0
          %             fprintf('(time %f) Progress %f\n', toc(tstart), counter/nchoosek(n,3));
          %           end
        end
      end
    end
    triangle_cols = repmat(1:3*nchoosek(n,3), 3, 1);
    vals = repmat([1; 1; -1], 3, nchoosek(n,3));
    triangle_cstr = sparse(triangle_idxs(:), triangle_cols(:), vals(:),...
      nchoosek(n,2), 3*nchoosek(n,3), 9*nchoosek(n,3))';
    
    % extract cost function from adjacency matrix
    E = zeros(nchoosek(n,2), 1);
    for i = 1:n
      for j = i+1:n
        E(coord_map(i, j)) = A(i, j);
      end
    end
    
    % Set up the LP
    fprintf('(time %f) Solving the LP.\n', toc(tstart));
    X = linprog(E, ... % cost
      -triangle_cstr, sparse([],[],[], size(triangle_cstr, 1), 1),... % inequality
      ones(1, nchoosek(n, 2)), 1, ... % equality
      zeros(nchoosek(n, 2), 1), inf * ones(nchoosek(n, 2), 1)); %nonnegativity
    %     fprintf('(time %f)\n', toc(tstart));
    
    fprintf('(time %f) OPT bound: %f.\n', toc(tstart), E'*X);
    
    % ROUNDING
    % Use Bourgain's embedding to get an \ell1 embedding with desired
    % distortion and not-too-large dimension. Then turn this \ell1
    % embedding into a cut metric using inductive algorithm.
    %
    % @TODO
  case 'sdp'
    % SEMIDEFINITE PROGRAMMING SOLVER
    % Use CVX for this.
    % SIZES: Storing the constraint matrices requires 3*nchoosek(n,3)
    % space, since they are +/- 1 matrices some compression is possible but
    % this gives worst case 4 GB space to store for n=1000 graph. So cannot
    % solve past this point at best.
    %
    EIG_TOL = 1e-4; % magnitude below this -> eigenvalue treated as zero
    fprintf('(time %f) Using the SDP solver.\n', toc(tstart));
    %     fprintf('(time %f) Creating constraint matrices.\n', toc(tstart));
    n = size(A,1);
    % there are nchoosek(n, 2) decision variables in the problem
    %     triangle_idxs = cell(nchoosek(n,3), 3);
    %     counter = 1;
    fprintf('(time %f) Inputting constraints, then solving the SDP.\n', toc(tstart));
    
    cvx_begin sdp
    %       cvx_precision high
    variable X(n,n) semidefinite
    % Objective/constraints for the sphere-constrained problem (not valid?)
    %       maximize trace(A*X)
    %       subject to
    %       diag(X) == ones(n,1)
    %       trace((ones(n,n) - spdiags(ones(n,1), 0, n, n)) * X) == n*(n-1) - n % (n+2)(n-1)/2
    % incorrectly formulated (below)...?
    %       trace( (ones(n,n) - X)*(ones(n,n) - eye(n))) == 1
    % Objective/constraints for the standard formulation of the problem
    % The objective double counts, so make the constraint double count
    % too for consistency.
    % This problem is all linear, hence homogeneous (to scalings of the
    % equality constraint).
    %       minimize trace(A*diag(X)*ones(1,n)) + trace(A*ones(n,1)*diag(X)') - 2*trace(A*X)
    minimize trace(D*X) - trace(A*X)
    subject to
    n*trace(X) - trace(ones(n,n)*X) == 0.5*n^2
    % Triangle inequality constraints
    for i = 1:n
      %         for j = 1:n
      %           for k = 1:n
      %         X(i,i) <= 1
      for j = i+1:n
        for k = j+1:n
          %             M1 = sparse(...
          %               [i; j; i],...
          %               [j; k; k],...
          %               [1; 1; -1],...
          %               n, n);
          %             M2 = sparse(...
          %               [i; j; i],...
          %               [k; k; j],...
          %               [1; 1; -1],...
          %               n, n);
          %             M3 = sparse(...
          %               [i; i; j],...
          %               [k; j; k],...
          %               [1; 1; -1],...
          %               n, n);
          %             trace(M1 * X) <= 1
          %             trace(M2 * X) <= 1
          %             trace(M3 * X) <= 1
          
          % These are for the sphere-constrained problem
          %               X(i,j) + X(j,k) - X(i,k) <= 1
          %               X(i,k) + X(j,k) - X(i,j) <= 1
          %               X(i,k) + X(i,j) - X(j,k) <= 1
          
          % These are for the standard problem
          X(i,j) + X(j,k) - X(i,k) - X(j,j) <= 0
          X(i,k) + X(j,k) - X(i,j) - X(k,k) <= 0
          X(i,k) + X(i,j) - X(j,k) - X(i,i) <= 0
          
          
          % triangle_idxs{counter, 1} = M1;
          % triangle_idxs{counter, 2} = M2;
          % triangle_idxs{counter, 3} = M3;
          %             counter = counter + 1;
          %           keyboard
        end
      end
    end
    %       for iter = 1:size(triangle_cstr,1)
    %         trace(triangle_cstr{iter} * X) <= 1
    %       end
    cvx_end
    fprintf('(time %f) Extracting metric from SDP solution.\n', toc(tstart));
    [V, E] = eig(X);
    eigens = diag(E);
    eigens(abs(eigens) <= EIG_TOL) = 0;
    metric_mx = V(:, eigens ~= 0) * diag(sqrt(eigens(eigens ~= 0)));
    distance_mx = pdistmx(metric_mx').^2; %\ell_2^2 distance matrix
    fprintf('(time %f) OPT bound: %f.\n', toc(tstart), trace(A*distance_mx)/n^2);
    % can also just take 2*ones(n,n) - 2*X
    
    
    
    
  case 'bm'
    %% Sloppy Burer-Monteiro-style solution of the problem
    % the ALM method subproblem seems to not have a nice closed-form
    % solution (if it has one it is not obvious), so this is going to run
    % extremely slow no matter what
    % here just solve the subproblems with a gradient descent solver
    % because it is hard to write out the system of linear equations we
    % need to solve... (and it requires O(n^4) space anyways!)
    
    % Step 1: prepare primitives
    % triangle inequality constraint matrices
    n = size(A,1);
    fprintf('(time %f) Using the BM solver.\n', toc(tstart));
    triangle_idxs = cell(nchoosek(n,3), 3);
    counter = 1;
    for i = 1:n
      for j = i+1:n
        for k = j+1:n
          M1 = sparse(...
            [i; j; i],...
            [j; k; k],...
            [1; 1; -1],...
            n, n);
          M2 = sparse(...
            [i; j; i],...
            [k; k; j],...
            [1; 1; -1],...
            n, n);
          M3 = sparse(...
            [i; i; j],...
            [k; j; k],...
            [1; 1; -1],...
            n, n);
          
          triangle_idxs{counter, 1} = M1 + M1';
          triangle_idxs{counter, 2} = M2 + M2';
          triangle_idxs{counter, 3} = M3 + M3';
          counter = counter + 1;
        end
      end
    end
    triangle_cstr = triangle_idxs(:);
    
    % laplacian
    L = D - A;
    % complete graph laplacian
    Lc = sparse(n*eye(n) - ones(n,n));
    
    % PARAMETERS
    lambda = 1; % ALM parameter
    p = 4;    % BM parameter
    alpha0 = 1;
    beta = 0.5;
    sigma = 1e-2;
    MAX_OUTER_ITER = 1e3;
    MAX_INNER_ITER = 1e2;
    INNER_TOL = 1e-4; % gradient norm tol
    OUTER_TOL = 1e-4; % objective decrease tol
    INFEAS_TOL = 1e-4;
    MIN_ALPHA = 0;
    % variables
    y0 = rand(nchoosek(n,3)*3, 1);
    kappa0 = rand(1, 1);
    mu0 = rand(nchoosek(n,3)*3, 1);
    Y0 = randn(n,p); Y0 = Y0 / sqrt(trace(Lc * (Y0*Y0'))); % initialize feasible
    
    mu = mu0;
    kappa = kappa0;
    y = y0;
    Y = Y0;
    % SOLVE PROBLEM via ALM METHOD
    fprintf('(time %f) Starting ALM.\n', toc(tstart));
    ofval = nan(MAX_OUTER_ITER, 1);
    infeas = nan(MAX_OUTER_ITER, 1);
    for idx_outer = 1:MAX_OUTER_ITER
      % compute the next (Y, y) value by finding minimizer of the
      % augmented lagrangian
      % use gradient.
      % initialize parameters.
      fval = nan(MAX_INNER_ITER, 1);
      % We start at the solution of the previous iteration
      alpha = alpha0;
      yprev = y;
      Yprev = Y;
      for idx_inner = 1:MAX_INNER_ITER
        X = Y + (idx_inner/(idx_inner+3))*(Y - Yprev);
%         x = y + (idx_inner/(idx_inner+3))*(y - yprev);
        x = y;
        
        accum1 = 0;
        accum2 = 0;
        XXt = X*X';
        
        for k = 1:size(triangle_cstr, 1)
          val = triangle_cstr{k} * X;
          accum1 = accum1 + mu(k) * val;
          accum2 = accum2 + (sum( vec(XXt .* triangle_cstr{k}) ) + x(k)) * val;
        end
        
        LcX = Lc*X;
        LX = L*X;
        % Gradient computations
        grad_X = 2 * (LX + kappa * LcX + accum1 ...
          + lambda * ( (sum( vec(Lc .* XXt)) - 1)*LcX + accum2));
        grad_x = zeros(size(mu,1), 1);
        for k = 1:length(grad_x)
          grad_x(k) = mu(k) + lambda * (sum( vec(XXt .* triangle_cstr{k}) ) + x(k));
        end
        
        % Gradient step
        % helpers for function evaluations
%         accum = @(x1)( sum( vec( x1 .* YYt) ) );
        f = @(x1, x2)( sum(vec(L .* (x1*x1'))) + kappa * ( sum(vec(Lc .* (x1*x1'))) - 1) ...
          + sum(x2 .* mu) + sum(cell_helper(triangle_cstr, x1*x1'))...
          + lambda/2 * (sum(vec(Lc .* (x1*x1'))) - 1)^2 ...
          + sum_square(cell_helper(triangle_cstr, x1*x1') + x2));
%         keyboard
        current_fval = f(X, x);
        while 1
          Z = X - alpha * grad_X;
          z = max(x - alpha * grad_x, 0); % projected gradient
          
          candidate_fval = f(Z, z);
          % backtracking line search
          if alpha <= MIN_ALPHA || current_fval - candidate_fval >= sigma*alpha*(norm(grad_X, 'fro')^2 + norm(grad_x, 2)^2)
            break;
          end
          
          alpha = alpha * beta;
%           fprintf('step size: %1.8f\n', alpha);
        end
        
        Yprev = Y;
        yprev = y;
        Y = Z;
        y = z;
        
        fval(idx_inner) = candidate_fval;
        if idx_inner > 1 && norm(grad_X, 'fro')^2 + norm(grad_x, 2)^2 <= INNER_TOL
          break;
        end
      end
%       keyboard
      % Apply the proximal operator on y
%       y = max(y, 0);
      fval(isnan(fval)) = [];
      fprintf('(time %f) completed inner iteration %d\n', toc(tstart), idx_outer);
%       keyboard
      
      % Outer iteration: update dual variables
      kappa = kappa + lambda * ( sum(vec(Lc .* (Y*Y'))) - 1);
      mu = mu + lambda * (y + cell_helper(triangle_cstr, Y*Y'));
      
%       f = @(x1, x2)( sum(vec(L .* (x1*x1'))) + kappa * ( sum(vec(Lc .* (x1*x1'))) - 1) ...
%         + sum(x2 .* mu) + sum(cell_helper(triangle_cstr, x1*x1'))...
%         + lambda/2 * (sum(vec(Lc .* (x1*x1'))) - 1)^2 ...
%         + sum_square(cell_helper(triangle_cstr, x1*x1') + x2));
%       ofval(idx_outer) = f(Y, y);
      ofval(idx_outer) = trace(L * (Y*Y'));
      infeas(idx_outer) = abs(trace(Lc*(Y*Y')) - 1)...
        + sum(abs(cell_helper(triangle_cstr, Y*Y') + y));
      fprintf('(time %f) ofval %f\n', toc(tstart), ofval(idx_outer));
      fprintf('(time %f) infeas %f\n', toc(tstart), infeas(idx_outer));
      
      if idx_outer > 1 && infeas(idx_outer - 1) - infeas(idx_outer) <= INFEAS_TOL
        % increase lambda if we didn't decrease enough in infeasibility
        lambda = lambda*2;
        fprintf('(time %f) increasing infeasibility parameter to %f\n', toc(tstart), lambda);
      end
      
      if idx_outer > 1 && abs(ofval(idx_outer - 1) - ofval(idx_outer)) <= OUTER_TOL ...
          && ofval(idx_outer - 1) - ofval(idx_outer) >= 0 && abs(infeas(idx_outer) - infeas(idx_outer-1)) <= INFEAS_TOL
        % quit if we decreased _and_ we decreased by a too-small amount.
        % and we're almost feasible
        break;
      end
    end
    
    keyboard
  otherwise
    error('unrecognized semimetric_mode setting');
end