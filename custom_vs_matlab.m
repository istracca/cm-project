% This script compares the performance of our MINRES and MATLAB's MINRES
% on graphs of different sizes, defined by exp_N_values (which represents
% the logarithm base 2 of the number of nodes in the graph) and
% E_N_ratio_values (which represents the ratio of the number of edges to
% the number of nodes in the graph). 
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


% Initialize LaTeX table string
table ='';

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
    % Custom MINRES:
    custom_times = zeros(1, number_of_instance);
    custom_relative_resvecs = zeros(1, number_of_instance);
    custom_iters = zeros(1, number_of_instance);

    % MATLAB MINRES:
    matlab_times = zeros(1, number_of_instance);
    matlab_relative_resvecs = zeros(1, number_of_instance);
    matlab_iters = zeros(1, number_of_instance);

    % Pre-allocate array to store condition numbers of A for different instances
    cond_nums = zeros(1, number_of_instance);

    % Cycle through different instances of same exp_N, E_N_ratio
    % configuration
    for instance = 1:number_of_instance
        % Create A from the dimacs file
        file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
        [D, E, y] = create_matrix_from_dimacs(file_path, flag, 1);
        A = sparse([D E'; E zeros(size(E,1))]);
        
        % Calculate the condition number of matrix A
        cond_nums(instance) = condest(A);

        % Custom MINRES
        % Start timer
        tic;
        [custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
        % Get the relative_resvec from the resvec
        custom_relative_resvec = custom_resvec / norm(y);
        % Stop timer
        custom_times(instance) = toc;
        % Final relative residual
        custom_relative_resvecs(instance) = custom_relative_resvec(end);
        % Adjust iteration count
        custom_iters(instance) = custom_iter-1;

        % MATLAB MINRES
        % Start timer
        tic;
        [x, flag, relres, iter, resvec] = minres(A, y, 1e-6, size(y,1));
        % Get the relative_resvec from the resvec
        relative_resvec = resvec / norm(y);
        % Remove the first entry for plotting
        relative_resvec = relative_resvec(2:end);
        % Stop timer
        matlab_times(instance) = toc;
        % Final relative residual
        matlab_relative_resvecs(instance) = relative_resvec(end);
        % Iteration count
        matlab_iters(instance) = iter;
    end

    % Calculate averages for custom MINRES
    % Average time
    mean_custom_time = mean(custom_times);
    % Average final relative residual
    mean_custom_resvec = mean(custom_relative_resvecs);
    % Average iterations
    mean_custom_iter = mean(custom_iters);
    
    % Calculate averages for MATLAB MINRES
    % Average time
    mean_matlab_time = mean(matlab_times);
    % Average final relative residual
    mean_matlab_resvec = mean(matlab_relative_resvecs);
    % Average iterations
    mean_matlab_iter = mean(matlab_iters);


    % Print results for this combination of parameters
    fprintf('Combination (%d, %d)\n', n, m);
    fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_custom_time, mean_custom_resvec, round(mean_custom_iter));
    fprintf('MATLAB MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_matlab_time, mean_matlab_resvec, round(mean_matlab_iter));
    fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_nums);
    % Append results to LaTeX table
    table = [table, sprintf('%d & %d & %d & %.4g & %d & %.4g & %.2g \\\\\n', n, m, round(mean_custom_iter), mean_custom_time, round(mean_matlab_iter), mean_matlab_time, mean_custom_time/mean_matlab_time)];
end

% Print the LaTeX table
fprintf('Table\n%s\n', table);
