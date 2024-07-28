% Plot of comparison between original, preconditioned and 
% sparse_preconditioned versions. A combination of exp_N, E_N_ratio
% parameters must be specified. For more information about these
% parameters, look at the main file "preconditioned_vs_original.m"

% Exponent to determine the number of nodes
exp_N = 8;
% Ratio to determine the size of the matrix m
E_N_ratio = 8;
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
cond_num = condest(A);
preconditioned_cond_num = condest((R')^(-1)*A*R^(-1));
sparse_preconditioned_cond_num = condest((sparse_R')^(-1)*A*sparse_R^(-1));

% Custom MINRES
% Start timer
tic;
% Run the algorithm
[custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
% Get the relative_resvec from the resvec
custom_relative_resvec = custom_resvec / norm(y);
% Stop timer
custom_time = toc;
% Final relative residual
custom_rel_resvec = custom_relative_resvec(end);
% Adjust iteration count
custom_iter = custom_iter-1;

% Custom MINRES preconditioned
% Start timer
tic;
% Run the algorithm
[preconditioned_x, preconditioned_flag, preconditioned_relres, preconditioned_iter, preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), R);
% Get the relative_resvec from the resvec
preconditioned_relative_resvec = preconditioned_resvec / norm(y);
% Stop timer
preconditioned_time = toc;
% Final relative residual
preconditioned_rel_resvec = preconditioned_relative_resvec(end);
% Adjust iteration count
preconditioned_iter = preconditioned_iter-1;

% Custom MINRES sparse preconditioned
% Start timer
tic;
% Run the algorithm
[sparse_preconditioned_x, sparse_preconditioned_flag, sparse_preconditioned_relres, sparse_preconditioned_iter, sparse_preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), sparse_R);
% Get the relative_resvec from the resvec
sparse_preconditioned_relative_resvec = sparse_preconditioned_resvec / norm(y);
% Stop timer
sparse_preconditioned_time = toc;
% Final relative residual
sparse_preconditioned_rel_resvec = sparse_preconditioned_relative_resvec(end);
% Adjust iteration count
sparse_preconditioned_iter = sparse_preconditioned_iter-1;


% Print results
fprintf('Combination (%d, %d)\n', n, m);
fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', custom_time, custom_rel_resvec, custom_iter);
fprintf('Condition number: %.4g\n', cond_num);
fprintf('Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', preconditioned_time, preconditioned_rel_resvec, preconditioned_iter);
fprintf('Preconditioned condition number: %.4g\n', preconditioned_cond_num);
fprintf('Sparse Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', sparse_preconditioned_time, sparse_preconditioned_rel_resvec, sparse_preconditioned_iter);
fprintf('Sparse Preconditioned condition number: %.4g\n', sparse_preconditioned_cond_num);



% Plot the relative_resvec for the three different cases
figure;
semilogy(custom_relative_resvec, 'r', 'LineWidth', 1.5);
hold on;
semilogy(preconditioned_relative_resvec, 'b', 'LineWidth', 1.5);
hold on;
semilogy(sparse_preconditioned_relative_resvec, 'g', 'LineWidth', 1.5);
hold off;

xlabel('Iterations');
ylabel('Relative Residual');
legend('Custom MINRES', 'Preconditioned custom MINRES', 'Sparse Preconditioned custom MINRES');
title(sprintf('Evolution of relative residuals for N=%d, E=%d', 2^exp_N, E_N_ratio*2^exp_N));
grid on;


% Plot the abs of positive and negative eigenvalues separately on a
% semilogarithmic scale, for the three different cases.

% Original system:
% Computation of eigenvalues
eigenvalues = real(eigs(A,n+m-1));

% Divide the eigenvalues into positive and negative
positive_eigenvalues = eigenvalues(eigenvalues > 0);
negative_eigenvalues = eigenvalues(eigenvalues < 0);

% Sort positive and negative eigenvalues separately
positive_sorted_eigenvalues = sort(positive_eigenvalues);
negative_sorted_eigenvalues = sort(-negative_eigenvalues);

figure;
% Plot of positive eigenvalues
subplot(2, 1, 1);
semilogy(positive_sorted_eigenvalues, 'b-', 'LineWidth', 2);
xlim([1,length(positive_sorted_eigenvalues)]);
ylabel('Positive eigenvalues');
title('Distribution of Positive Eigenvalues');
grid on;

% Plot of abs of negative eigenvalues
subplot(2, 1, 2);
semilogy(negative_sorted_eigenvalues, 'b-', 'LineWidth', 2); % Plotta gli autovalori come linea continua bluxlabel('Eigenvalues (Positive)');
xlim([1,length(negative_sorted_eigenvalues)]);
ylabel('Abs(Negative eigenvalues)');
title('Distribution of Negative Eigenvalues');
grid on;

sgtitle('Eigenvalues of A');


% Preconditioned version:
% Computation of eigenvalues
eigenvalues = real(eigs(R'^(-1)*A*R^(-1),n+m-1));

% Divide the eigenvalues into positive and negative
positive_eigenvalues = eigenvalues(eigenvalues > 0);
negative_eigenvalues = eigenvalues(eigenvalues < 0);

% Sort positive and negative eigenvalues separately
positive_sorted_eigenvalues = sort(positive_eigenvalues);
negative_sorted_eigenvalues = sort(-negative_eigenvalues);

figure;
% Plot of positive eigenvalues
subplot(2, 1, 1);
semilogy(positive_sorted_eigenvalues, 'b-', 'LineWidth', 2);
xlim([1,length(positive_sorted_eigenvalues)]);
ylabel('Positive eigenvalues');
title('Distribution of Positive Eigenvalues');
grid on;

% Plot of abs of negative eigenvalues
subplot(2, 1, 2);
semilogy(negative_sorted_eigenvalues, 'b-', 'LineWidth', 2);
xlim([1,length(negative_sorted_eigenvalues)]);
ylabel('Abs(Negative eigenvalues)');
title('Distribution of Negative Eigenvalues');
grid on;

sgtitle('Eigenvalues in Preconditioned Case');


% Sparse preconditioned version:
% Computation of eigenvalues
eigenvalues = real(eigs((sparse_R^(-1))'*A*sparse_R^(-1),n+m-1));

% Divide the eigenvalues into positive and negative
positive_eigenvalues = eigenvalues(eigenvalues > 0);
negative_eigenvalues = eigenvalues(eigenvalues < 0);

% Sort positive and negative eigenvalues separately
positive_sorted_eigenvalues = sort(positive_eigenvalues);
negative_sorted_eigenvalues = sort(-negative_eigenvalues);

figure;
% Plot of positive eigenvalues
subplot(2, 1, 1);
semilogy(positive_sorted_eigenvalues, 'b-', 'LineWidth', 2);
xlim([1,length(positive_sorted_eigenvalues)]);
ylabel('Positive eigenvalues');
title('Distribution of Positive Eigenvalues');
grid on;

% Plot of abs of negative eigenvalues
subplot(2, 1, 2);
semilogy(negative_sorted_eigenvalues, 'b-', 'LineWidth', 2);
xlim([1,length(negative_sorted_eigenvalues)]);
ylabel('Abs(Negative eigenvalues)');
title('Distribution of Negative Eigenvalues');
grid on;

sgtitle('Eigenvalues in Sparse Preconditioned Case');