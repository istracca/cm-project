% The script reads a matrix instance from a file, constructs a sparse matrix A, 
% computes the condition number of A, and then applies both custom_minres and MATLAB's minres 
% to solve the linear system A*x = y. Finally, it creates a plot comparing relative residuals 
% at each iteration for both methods.
% Outputs:
% - Printed results for custom_minres and MATLAB's minres (time, residual, iterations)
% - Plot comparing the evolution of relative residuals for both methods

% Exponent to determine the number of nodes
exp_N = 8;
% Ratio to determine the size of the matrix m
E_N_ratio = 16;
% Number of instances to average the results over for each combination of
% parameters
instance = 3;
% Define the flag for generation of D
flag = 3;
% You can change the value of the variables above to obtain the results for
% different parameters' combinations.


% Number of nodes
n = 2^exp_N;
% Number of edges
m = n*E_N_ratio;

% Create A from the dimacs file
file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
[D, E, y] = create_matrix_from_dimacs(file_path, flag, 1);
A = sparse([D E'; E zeros(size(E,1))]);

% Calculate the condition number of matrix A
cond_num = condest(A);

% Custom MINRES
% Start timer
tic;
% Run custom MINRES
[custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
% Get the relative_resvec from the resvec
custom_relative_resvec = custom_resvec / norm(y);
% Stop timer
custom_time = toc;
% Final relative residual
custom_resvec = custom_relative_resvec(end);
% Adjust iteration count
custom_iter = custom_iter-1;

% MATLAB MINRES
% Start timer
tic;
% Run MATLAB MINRES
[x, flag, relres, iter, resvec] = minres(A, y, 1e-6, size(y,1));
% Normalize residuals
relative_resvec = resvec / norm(y);
% Remove the first entry for plotting
relative_resvec = relative_resvec(2:end);
% Stop timer
matlab_time = toc;
% Final relative residual
matlab_resvec = relative_resvec(end);
% Iteration count
matlab_iter = iter;

% Print results
fprintf('Parameters Combination (%d, %d)\n', n, m);
fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', custom_time, custom_resvec, custom_iter);
fprintf('MATLAB MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', matlab_time, matlab_resvec, matlab_iter);
fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_num);

% Plot the evolution of relative residuals for Custom and MATLAB MINRES
figure;
semilogy(custom_relative_resvec, 'b', 'LineWidth', 1.5);
hold on;
semilogy(relative_resvec, 'r', 'LineWidth', 1.5);
hold off;

xlabel('Iterations');
ylabel('Relative Residual');
legend('Custom MINRES', 'Matlab MINRES');
title(sprintf('Evolution of relative residuals for N=%d, E=%d', 2^exp_N, E_N_ratio*2^exp_N));
grid on;
