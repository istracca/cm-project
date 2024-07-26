% Plot di confronto tra versione originale, preconditioned e
% sparse_preconditioned. Bisogna specificare una combinazione di parametri
% del grafo

exp_N = 8;
E_N_ratio = 16;
instance = 1;
n = 2^exp_N;
m = n*E_N_ratio;

file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
[D, E, y] = create_matrix_from_dimacs(file_path, 1);
A = sparse([D E'; E zeros(size(E,1))]);

schur_component = E*D^(-1)*E';
[L_schur, M_schur] = ldl(schur_component);
L = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) L_schur]);
M = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) M_schur]);
R = sparse(M*L');

schur_component_trimmed = trim_schur_component(schur_component,1e-5,0);
[sparse_L_schur, sparse_M_schur] = ldl(schur_component_trimmed);
sparse_L = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) sparse_L_schur]);
sparse_M = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) sparse_M_schur]);
sparse_R = sparse(sparse_M*sparse_L');

% Calcola il numero di condizionamento della matrice A
cond_num = condest(A);
preconditioned_cond_num = condest((R')^(-1)*A*R^(-1));
sparse_preconditioned_cond_num = condest((sparse_R')^(-1)*A*sparse_R^(-1));

% Custom MINRES
tic;
[custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
custom_relative_resvec = custom_resvec / norm(y);
custom_time = toc;
custom_resvec = custom_relative_resvec(end);
custom_iter = custom_iter-1;

% Custom MINRES preconditioned
tic;
[preconditioned_x, preconditioned_flag, preconditioned_relres, preconditioned_iter, preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), R);
preconditioned_relative_resvec = preconditioned_resvec / norm(y);
preconditioned_time = toc;
preconditioned_resvec = preconditioned_relative_resvec(end);
preconditioned_iter = preconditioned_iter-1;

% Custom MINRES sparse preconditioned
tic;
[sparse_preconditioned_x, sparse_preconditioned_flag, sparse_preconditioned_relres, sparse_preconditioned_iter, sparse_preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), sparse_R);
sparse_preconditioned_relative_resvec = sparse_preconditioned_resvec / norm(y);
sparse_preconditioned_time = toc;
sparse_preconditioned_resvec = sparse_preconditioned_relative_resvec(end);
sparse_preconditioned_iter = sparse_preconditioned_iter-1;


fprintf('Combinazione (%d, %d)\n', n, m);
fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', custom_time, custom_resvec, custom_iter);
fprintf('Condition number: %.4g\n', cond_num);
fprintf('Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', preconditioned_time, preconditioned_resvec, preconditioned_iter);
fprintf('Preconditioned condition number: %.4g\n', preconditioned_cond_num);
fprintf('Sparse Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', sparse_preconditioned_time, sparse_preconditioned_resvec, sparse_preconditioned_iter);
fprintf('Sparse Preconditioned condition number: %.4g\n', sparse_preconditioned_cond_num);




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

% Calcolo degli autovalori
eigenvalues = real(eigs(A,2303));

% Separiamo gli autovalori in positivi e negativi
positive_eigenvalues = eigenvalues(eigenvalues > 0);
negative_eigenvalues = eigenvalues(eigenvalues < 0);

positive_sorted_eigenvalues = sort(positive_eigenvalues);
negative_sorted_eigenvalues = sort(-negative_eigenvalues);

% Crea una figura per il grafico degli autovalori positivi
figure;

% Plotta la stima della densità di kernel per gli autovalori positivi
subplot(2, 1, 1);
semilogy(positive_sorted_eigenvalues, 'b-', 'LineWidth', 2); % Plotta gli autovalori come linea continua bluxlabel('Eigenvalues (Positive)');
ylabel('Positive eigenvalues');
title('Distribution of Positive Eigenvalues');
set(gca, 'XScale', 'log'); % Imposta l'asse x su scala logaritmica
grid on;

% Crea una figura per il grafico degli autovalori negativi
subplot(2, 1, 2);
semilogy(negative_sorted_eigenvalues, 'b-', 'LineWidth', 2); % Plotta gli autovalori come linea continua bluxlabel('Eigenvalues (Positive)');
ylabel('Abs(Negative eigenvalues)');
title('Distribution of Negative Eigenvalues');
set(gca, 'XScale', 'log'); % Imposta l'asse x su scala logaritmica
grid on;

% Aggiunge un titolo generale alla figura
sgtitle('Eigenvalues of A');