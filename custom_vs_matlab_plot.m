% Plot di confronto tra l'evoluzione dei relative_residuals per 
% matlab_minres e custom_minres. Bisogna specificare una combinazione di
% parametri del grafo

exp_N = 8;
E_N_ratio = 16;
instance = 3;
n = 2^exp_N;
m = n*E_N_ratio;

file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
[D, E, y] = create_matrix_from_dimacs(file_path, 1);
A = sparse([D E'; E zeros(size(E,1))]);

% Calcola il numero di condizionamento della matrice A
cond_num = condest(A);

% Custom MINRES
tic;
[custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
custom_relative_resvec = custom_resvec / norm(y);
custom_time = toc;
custom_resvec = custom_relative_resvec(end);
custom_iter = custom_iter-1;

% MATLAB MINRES
tic;
[x, flag, relres, iter, resvec] = minres(A, y, 1e-6, size(y,1));
relative_resvec = resvec / norm(y);
relative_resvec = relative_resvec(2:end);
matlab_time = toc;
matlab_resvec = relative_resvec(end);
matlab_iter = iter;

% Stampa i risultati per questa combinazione
fprintf('Combinazione (%d, %d)\n', n, m);
fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', custom_time, custom_resvec, custom_iter);
fprintf('MATLAB MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', matlab_time, matlab_resvec, matlab_iter);
fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_num);

% Plot dell'evoluzione dei residui relativi per l'ultima combinazione
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
