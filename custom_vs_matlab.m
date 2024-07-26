% In questo file vengono confrontati matlab_minres e custom_minres, su
% grafici di dimensioni definite da exp_N_values (che rappresenta il
% logaritmo in base 2 del numero di nodi del grafo) e E_N_ratio_values (che
% rappresenta il rapporto tra numero di archi e numero di nodi del grafo)
% number_of_instance specifica invece su quante istanze di grafi, per ogni
% combinazione, mediare i risultati ottenuti.
% Viene stampata anche la tabella con la sintassi di LateX, per nostra
% comodità

exp_N_values = [8, 8, 8, 10, 10, 10, 12, 12];
E_N_ratio_values = [8, 16, 32, 8, 32, 64, 8, 64];
number_of_instance = 3;

table ='';

for idx = 1:length(exp_N_values)
    exp_N = exp_N_values(idx);
    E_N_ratio = E_N_ratio_values(idx);
    n = 2^exp_N;
    m = n*E_N_ratio;

    custom_times = zeros(1, number_of_instance);
    custom_resvecs = zeros(1, number_of_instance);
    custom_iters = zeros(1, number_of_instance);

    matlab_times = zeros(1, number_of_instance);
    matlab_resvecs = zeros(1, number_of_instance);
    matlab_iters = zeros(1, number_of_instance);

    cond_nums = zeros(1, number_of_instance);

    for instance = 1:number_of_instance
        file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
        [D, E, y] = create_matrix_from_dimacs(file_path, 1);
        A = sparse([D E'; E zeros(size(E,1))]);
        
        % Calcola il numero di condizionamento della matrice A
        cond_nums(instance) = condest(A);

        % Custom MINRES
        tic;
        [custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, size(y,1));
        custom_relative_resvec = custom_resvec / norm(y);
        custom_times(instance) = toc;
        custom_resvecs(instance) = custom_relative_resvec(end);
        custom_iters(instance) = custom_iter-1;

        % MATLAB MINRES
        tic;
        [x, flag, relres, iter, resvec] = minres(A, y, 1e-6, size(y,1));
        relative_resvec = resvec / norm(y);
        relative_resvec = relative_resvec(2:end);
        matlab_times(instance) = toc;
        matlab_resvecs(instance) = relative_resvec(end);
        matlab_iters(instance) = iter;
    end

    % Calcola le medie
    mean_custom_time = mean(custom_times);
    mean_custom_resvec = mean(custom_resvecs);
    mean_custom_iter = mean(custom_iters);

    mean_matlab_time = mean(matlab_times);
    mean_matlab_resvec = mean(matlab_resvecs);
    mean_matlab_iter = mean(matlab_iters);


    % Stampa i risultati per questa combinazione
    fprintf('Combinazione (%d, %d)\n', n, m);
    fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_custom_time, mean_custom_resvec, round(mean_custom_iter));
    fprintf('MATLAB MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_matlab_time, mean_matlab_resvec, round(mean_matlab_iter));
    fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_nums);
    table = [table, sprintf('%d & %d & %d & %.4g & %d & %.4g & %.2g \\\\\n', n, m, round(mean_custom_iter), mean_custom_time, round(mean_matlab_iter), mean_matlab_time, mean_custom_time/mean_matlab_time)];
end

fprintf('Tabella\n%s\n', table);

% 
% % Plot dell'evoluzione dei residui relativi per l'ultima combinazione
% figure;
% semilogy(custom_relative_resvec, 'b', 'LineWidth', 1.5);
% hold on;
% semilogy(relative_resvec, 'r', 'LineWidth', 1.5);
% hold off;
% 
% xlabel('Iterazioni');
% ylabel('Relative Residual');
% legend('Custom MINRES', 'Matlab MINRES');
% title(sprintf('Evolution of relative residuals for N=%d, E=%d', 2^exp_N, E_N_ratio*2^exp_N));
% grid on;
