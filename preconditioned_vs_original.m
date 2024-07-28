% This file compares custom_minres, version with preconditioner and version
% with sparse approximation of preconditioner, on graphs of sizes defined 
% by exp_N_values (representing the logarithm to base 2 of the number of 
% nodes in the graph) and E_N_ratio_values (representing the ratio of the 
% number of edges to the number of nodes in the graph).
% number_of_instance, on the other hand, specifies over how many graph 
% instances, for each combination, the results obtained should be averaged.
% The script also prints a table in LaTeX syntax with the results.

% Vector of exponents defining the size of the matrices (n = 2^exp_N)
exp_N_values = [8, 8, 8, 10, 10, 10, 12, 12];
% Vector of ratios defining the number of edges (m = n * E_N_ratio)
E_N_ratio_values = [8, 16, 32, 8, 32, 64, 8, 64];
% Number of instances to average the results over for each combination of
% parameters
number_of_instance = 3;
% Define the flag for generation of D
flag = 3;

% Initialize the two LaTeX strings
first_table = '';
second_table = '';


% Cycle through all the combinations of exp_N_values, E_N_ratio_values
for idx = 1:length(exp_N_values)
    % Current exponent value
    exp_N = exp_N_values(idx);
    % Current edge-to-node ratio
    E_N_ratio = E_N_ratio_values(idx);
    % Number of nodes
    n = 2^exp_N;
    % Number of edges
    m = n*E_N_ratio;

    % Pre-allocate the arrays to store the times of execution, the resvecs
    % and the number of iterations for different instances.
    % Original system:
    custom_times = zeros(1, number_of_instance);
    custom_relative_resvecs = zeros(1, number_of_instance);
    custom_iters = zeros(1, number_of_instance);

    % Preconditioned version:
    preconditioned_times = zeros(1, number_of_instance);
    preconditioned_relative_resvecs = zeros(1, number_of_instance);
    preconditioned_iters = zeros(1, number_of_instance);

    % Sparse preconditioned version:
    sparse_preconditioned_times = zeros(1, number_of_instance);
    sparse_preconditioned_relative_resvecs = zeros(1, number_of_instance);
    sparse_preconditioned_iters = zeros(1, number_of_instance);

    % Pre-allocate the arrays to store the condition numbers in the three
    % cases, for different instances.
    cond_nums = zeros(1, number_of_instance);
    preconditioned_cond_nums = zeros(1, number_of_instance);
    sparse_preconditioned_cond_nums = zeros(1, number_of_instance);

    % Cycle through different instances of same exp_N, E_N_ratio
    % configuration
    for instance = 1:number_of_instance
        % Create A from the dimacs file
        file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
        [D, E, y] = create_matrix_from_dimacs(file_path, flag, 1);
        A = sparse([D E'; E zeros(size(E,1))]);

        % Preconditioned version
        schur_component = E*D^(-1)*E';
        [L_schur, M_schur] = ldl(schur_component);
        L = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) L_schur]);
        M = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) M_schur]);
        R = sparse(M*L');

        % Sparse preconditioned version
        % k=0 in trim_schur_component implies that we are removing all 
        % non-diagonal entries of schur_component
        schur_component_trimmed = trim_schur_component(schur_component,1e-5,0);
        [sparse_L_schur, sparse_M_schur] = ldl(schur_component_trimmed);
        sparse_L = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) sparse_L_schur]);
        sparse_M = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) sparse_M_schur]);
        sparse_R = sparse(sparse_M*sparse_L');
        
        % Compute the condition numbers of the three version matrices
        cond_nums(instance) = condest(A);
        preconditioned_cond_nums(instance) = condest((R')^(-1)*A*R^(-1));
        sparse_preconditioned_cond_nums(instance) = condest((sparse_R')^(-1)*A*sparse_R^(-1));

        % Custom MINRES
        tic;
        % Run the algorithm
        [custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
        % Get the relative_resvec from the resvec
        custom_relative_resvec = custom_resvec / norm(y);
        % Stop timer
        custom_times(instance) = toc;
        % Final relative residual
        custom_relative_resvecs(instance) = custom_relative_resvec(end);
        % Adjust iteration count
        custom_iters(instance) = custom_iter-1;

        % Preconditioned version
        tic;
        % Run the algorithm
        [preconditioned_x, preconditioned_flag, preconditioned_relres, preconditioned_iter, preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), R);
        % Get the relative_resvec from the resvec
        preconditioned_relative_resvec = preconditioned_resvec / norm(y);
        % Stop timer
        preconditioned_times(instance) = toc;
        % Final relative residual
        preconditioned_relative_resvecs(instance) = preconditioned_relative_resvec(end);
        % Adjust iteration count
        preconditioned_iters(instance) = preconditioned_iter-1;

        % Sparse preconditioned version
        tic;
        % Run the algorithm
        [sparse_preconditioned_x, sparse_preconditioned_flag, sparse_preconditioned_relres, sparse_preconditioned_iter, sparse_preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), sparse_R);
        % Get the relative_resvec from the resvec
        sparse_preconditioned_relative_resvec = sparse_preconditioned_resvec / norm(y);
        % Stop timer
        sparse_preconditioned_times(instance) = toc;
        % Final relative residual
        sparse_preconditioned_relative_resvecs(instance) = sparse_preconditioned_relative_resvec(end);
        % Adjust iteration count
        sparse_preconditioned_iters(instance) = sparse_preconditioned_iter-1;

    end

    % Compute the averages
    mean_custom_time = mean(custom_times);
    mean_custom_resvec = mean(custom_relative_resvecs);
    mean_custom_iter = mean(custom_iters);
    mean_custom_cond_nums = mean(cond_nums);

    mean_preconditioned_time = mean(preconditioned_times);
    mean_preconditioned_resvec = mean(preconditioned_relative_resvecs);
    mean_preconditioned_iter = mean(preconditioned_iters);
    mean_preconditioned_cond_nums = mean(preconditioned_cond_nums);

    mean_sparse_preconditioned_time = mean(sparse_preconditioned_times);
    mean_sparse_preconditioned_resvec = mean(sparse_preconditioned_relative_resvecs);
    mean_sparse_preconditioned_iter = mean(sparse_preconditioned_iters);
    mean_sparse_preconditioned_cond_nums = mean(sparse_preconditioned_cond_nums);


    % Print results of the exp_N, E_N_ratio combination
    fprintf('Combination (%d, %d)\n', n, m);
    fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_custom_time, mean_custom_resvec, round(mean_custom_iter));
    fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_nums);
    fprintf('Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_preconditioned_time, mean_preconditioned_resvec, round(mean_preconditioned_iter));
    fprintf('Preconditioned condition numbers: %.4g, %.4g, %.4g\n', preconditioned_cond_nums);
    fprintf('Sparse Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_sparse_preconditioned_time, mean_sparse_preconditioned_resvec, round(mean_sparse_preconditioned_iter));
    fprintf('Sparse Preconditioned condition numbers: %.4g, %.4g, %.4g\n', sparse_preconditioned_cond_nums);
    % Update the two LaTeX tables
    first_table = [first_table, sprintf('%d & %d & %d & %.4g &  %.4g & %.4g & %.4g & %d & %.4g & %.4g & %.4g & %.4g \\\\\n', n, m, round(mean_custom_iter), mean_custom_resvec, cond_nums, round(mean_preconditioned_iter), mean_preconditioned_resvec, preconditioned_cond_nums)];
    second_table = [second_table, sprintf('%d & %d & %d & %.4g & %.4g & %d & %.4g & %.4g & %d & %.4g & %.4g \\\\\n', n, m, round(mean_custom_iter), mean_custom_time, mean_custom_cond_nums, round(mean_preconditioned_iter), mean_preconditioned_time, mean_preconditioned_cond_nums, round(mean_sparse_preconditioned_iter),mean_sparse_preconditioned_time,mean_sparse_preconditioned_cond_nums)];
end

% Print the two LaTeX tables
fprintf('Prima tabella\n%s\n', first_table);
fprintf('Seconda tabella\n%s\n', second_table);