exp_N_values = [8, 8, 8, 10, 10, 10];
E_N_ratio_values = [8, 16, 32, 8, 32, 64];
preconditioning_ratios = [0,0.01,0.05,0.1];
number_of_instance = 3;
warning('off', 'all');

% Array per salvare i risultati
results = zeros(length(exp_N_values), 9);
first_table = '';
second_table = '';
condition_numbers = zeros(length(exp_N_values), number_of_instance);
preconditioned_cond_numbers = zeros(length(exp_N_values), number_of_instance);
sparse_preconditioned_cond_numbers = zeros(length(exp_N_values), number_of_instance);

for idx = 1:length(exp_N_values)
    exp_N = exp_N_values(idx);
    E_N_ratio = E_N_ratio_values(idx);
    n = 2^exp_N;
    m = n*E_N_ratio;

    for j = 1:klength(preconditioning_ratios)
        times(j) = zeros(1, number_of_instance);
        resvecs(j) = zeros(1, number_of_instance);
        iters(j) = zeros(1, number_of_instance);
        cond_nums(j) = zeros(1,number_of_instance);

    for instance = 1:number_of_instance
        file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
        [D, E, y] = create_matrix_from_dimacs(file_path, 1);
        A = sparse([D E'; E zeros(size(E,1))]);

        schur_component = E*D^(-1)*E';
        [L_schur, M_schur] = ldl(schur_component);
        L = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) L_schur]);
        M = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) M_schur]);
        R = sparse(M*L');


        for j=1:length(preconditioning_ratios)
            schur_component_trimmed(j) = trim_schur_component(schur_component,1e-5,preconditioning_ratios(j)*(n-1)*(n-1));
            [sparse_L_schur(j), sparse_M_schur(j)] = ldl(schur_component_trimmed(j));
            sparse_L(j) = sparse([spdiags(ones(m,1),0,m,m) zeros(m,(n-1)); zeros((n-1),m) sparse_L_schur(j)]);
            sparse_M(j) = sparse([sqrt(D) zeros(m,(n-1)); zeros((n-1),m) sparse_M_schur(j)]);
            sparse_R(j) = sparse(sparse_M(j)*sparse_L(j)');
        
            % Calcola il numero di condizionamento della matrice A
            cond_nums(j)(instance) = condest(A);
            preconditioned_cond_nums(instance) = condest((R')^(-1)*A*R^(-1));
            sparse_preconditioned_cond_nums(instance) = condest((sparse_R')^(-1)*A*sparse_R^(-1));

        % Custom MINRES
        tic;
        [custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
        custom_relative_resvec = custom_resvec / norm(y);
        custom_times(instance) = toc;
        custom_resvecs(instance) = custom_relative_resvec(end);
        custom_iters(instance) = custom_iter-1;

        % Custom MINRES preconditioned
        tic;
        [preconditioned_x, preconditioned_flag, preconditioned_relres, preconditioned_iter, preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), R);
        preconditioned_relative_resvec = preconditioned_resvec / norm(y);
        preconditioned_times(instance) = toc;
        preconditioned_resvecs(instance) = preconditioned_relative_resvec(end);
        preconditioned_iters(instance) = preconditioned_iter-1;

        % Custom MINRES sparse preconditioned
        tic;
        [sparse_preconditioned_x, sparse_preconditioned_flag, sparse_preconditioned_relres, sparse_preconditioned_iter, sparse_preconditioned_resvec] = custom_minres_preconditioned(A, y, 1e-6, size(y,1), sparse_R);
        sparse_preconditioned_relative_resvec = sparse_preconditioned_resvec / norm(y);
        sparse_preconditioned_times(instance) = toc;
        sparse_preconditioned_resvecs(instance) = sparse_preconditioned_relative_resvec(end);
        sparse_preconditioned_iters(instance) = sparse_preconditioned_iter-1;

    end

    % Calcola le medie
    mean_custom_time = mean(custom_times);
    mean_custom_resvec = mean(custom_resvecs);
    mean_custom_iter = mean(custom_iters);
    mean_custom_cond_nums = mean(cond_nums);

    mean_preconditioned_time = mean(preconditioned_times);
    mean_preconditioned_resvec = mean(preconditioned_resvecs);
    mean_preconditioned_iter = mean(preconditioned_iters);
    mean_preconditioned_cond_nums = mean(preconditioned_cond_nums);

    mean_sparse_preconditioned_time = mean(sparse_preconditioned_times);
    mean_sparse_preconditioned_resvec = mean(sparse_preconditioned_resvecs);
    mean_sparse_preconditioned_iter = mean(sparse_preconditioned_iters);
    mean_sparse_preconditioned_cond_nums = mean(sparse_preconditioned_cond_nums);


    % Salva i risultati
    results(idx, :) = [mean_custom_time, mean_custom_resvec, mean_custom_iter, mean_preconditioned_time, mean_preconditioned_resvec, mean_preconditioned_iter, mean_sparse_preconditioned_time, mean_sparse_preconditioned_resvec, mean_sparse_preconditioned_iter];


    % Stampa i risultati per questa combinazione
    fprintf('Combinazione (%d, %d)\n', n, m);
    fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_custom_time, mean_custom_resvec, round(mean_custom_iter));
    fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_nums);
    fprintf('Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_preconditioned_time, mean_preconditioned_resvec, round(mean_preconditioned_iter));
    fprintf('Preconditioned condition numbers: %.4g, %.4g, %.4g\n', preconditioned_cond_nums);
    fprintf('Sparse Preconditioned MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_sparse_preconditioned_time, mean_sparse_preconditioned_resvec, round(mean_sparse_preconditioned_iter));
    fprintf('Sparse Preconditioned condition numbers: %.4g, %.4g, %.4g\n', sparse_preconditioned_cond_nums);
    first_table = [first_table, sprintf('%d & %d & %d & %.4g &  %.4g & %.4g & %.4g & %d & %.4g & %.4g & %.4g & %.4g \\\\\n', n, m, round(mean_custom_iter), mean_custom_resvec, cond_nums, round(mean_preconditioned_iter), mean_preconditioned_resvec, preconditioned_cond_nums)];
    second_table = [second_table, sprintf('%d & %d & %d & %.4g & %.4g & %d & %.4g & %.4g & %d & %.4g & %.4g \\\\\n', n, m, round(mean_custom_iter), mean_custom_time, mean_custom_cond_nums, round(mean_preconditioned_iter), mean_preconditioned_time, mean_preconditioned_cond_nums, round(mean_sparse_preconditioned_iter),mean_sparse_preconditioned_time,mean_sparse_preconditioned_cond_nums)];
end

fprintf('Prima tabella\n%s\n', first_table);
fprintf('Seconda tabella\n%s\n', second_table);

% Plot dell'evoluzione dei residui relativi per l'ultima combinazione
% figure;
% semilogy(custom_relative_resvec, 'r', 'LineWidth', 1.5);
% hold on;
% semilogy(preconditioned_relative_resvec, 'b', 'LineWidth', 1.5);
% hold on;
% semilogy(sparse_preconditioned_relative_resvec, 'g', 'LineWidth', 1.5);
% hold off;
% 
% xlabel('Iterazioni');
% ylabel('Relative Residual');
% legend('Custom MINRES', 'Preconditioned custom MINRES', 'Sparse Preconditioned custom MINRES');
% title(sprintf('Evolution of relative residuals for N=%d, E=%d', 2^exp_N, E_N_ratio*2^exp_N));
% grid on;

% eigenvalues=real(eigs((R')^(-1)*A*R^(-1),2303));
% 
% % Calcola la stima della densità di kernel
% [f, xi] = ksdensity(eigenvalues);

% % Crea una figura per il grafico
% figure;
% 
% % Plotta la stima della densità di kernel
% plot(xi, f, 'LineWidth', 2);
% 
% % Aggiungi etichette agli assi e un titolo
% xlabel('Autovalori');
% ylabel('Densità');
% title('Distribuzione degli Autovalori (Scala Logaritmica)');
% 
% % Imposta l'asse x su scala logaritmica
% set(gca, 'XScale', 'log');
% 
% % Aggiungi una griglia al grafico
% grid on;

% kurt = kurtosis(eigenvalues);
% 
% % Visualizza il risultato
% disp(['Kurtosis degli autovalori: ', num2str(kurt)]);
% 
% epsilon = 1e-10;  % Piccola costante positiva
% f = f + epsilon;
% 
% % Calcola l'indice di Shannon
% shannon_index = -sum(f .* log(f));
% 
% % Mostra l'indice di Shannon
% disp(['Indice di Shannon: ', num2str(shannon_index)]);