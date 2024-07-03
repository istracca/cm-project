exp_N_values = [8, 8, 8, 10, 10, 10, 12, 12];
E_N_ratio_values = [8, 16, 32, 8, 32, 64, 8, 64];
number_of_instance = 3;

% Array per salvare i risultati
results = zeros(length(exp_N_values), 6); % 6 columns for mean values of the results
condition_numbers = zeros(length(exp_N_values), 3); % 3 columns for condition numbers of each instance

for idx = 1:length(exp_N_values)
    exp_N = exp_N_values(idx);
    E_N_ratio = E_N_ratio_values(idx);

    custom_times = zeros(1, number_of_instance);
    custom_resvecs = zeros(1, number_of_instance);
    custom_iters = zeros(1, number_of_instance);

    matlab_times = zeros(1, number_of_instance);
    matlab_resvecs = zeros(1, number_of_instance);
    matlab_iters = zeros(1, number_of_instance);

    cond_nums = zeros(1, number_of_instance);

    for instance = 1:number_of_instance
        file_path = sprintf('graph_instances/net%d_%d_%d.dmx', exp_N, E_N_ratio, instance);
        [A, y, complement_eigenvalues, D_values] = create_matrix_from_dimacs(file_path, 1);
        
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

    % Salva i risultati
    results(idx, :) = [mean_custom_time, mean_custom_resvec, mean_custom_iter, mean_matlab_time, mean_matlab_resvec, mean_matlab_iter];

    % Salva i numeri di condizionamento
    condition_numbers(idx, :) = cond_nums;

    % Stampa i risultati per questa combinazione
    fprintf('Combinazione (%d, %d)\n', exp_N, E_N_ratio);
    fprintf('Custom MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_custom_time, mean_custom_resvec, round(mean_custom_iter));
    fprintf('MATLAB MINRES - Time: %.6f seconds, ResVec: %.6g, Iterations: %d\n', mean_matlab_time, mean_matlab_resvec, round(mean_matlab_iter));
    fprintf('Condition numbers: %.4g, %.4g, %.4g\n', cond_nums);
end

% Crea una tabella con i risultati
results_table = array2table(results, 'VariableNames', ...
    {'Custom_Time', 'Custom_ResVec', 'Custom_Iter', 'Matlab_Time', 'Matlab_ResVec', 'Matlab_Iter'});
results_table.Exp_N = exp_N_values';
results_table.E_N_Ratio = E_N_ratio_values';

% Crea una tabella con i numeri di condizionamento
condition_table = array2table(condition_numbers, 'VariableNames', ...
    {'Cond_Num_1', 'Cond_Num_2', 'Cond_Num_3'});
condition_table.Exp_N = exp_N_values';
condition_table.E_N_Ratio = E_N_ratio_values';

% Combina le due tabelle
final_table = [results_table, condition_table(:, 3:end)];

% Salva la tabella in un file CSV
writetable(final_table, 'minres_results_with_cond_numbers.csv');

% Plot dell'evoluzione dei residui relativi per l'ultima combinazione
figure;
semilogy(custom_relative_resvec, 'b', 'LineWidth', 1.5);
hold on;
semilogy(relative_resvec, 'r', 'LineWidth', 1.5);
hold off;

xlabel('Iterazioni');
ylabel('Relative Residual');
legend('Custom MINRES', 'Matlab MINRES');
title(sprintf('Evolution of relative residuals for N=%d, E=%d', 2^exp_N, E_N_ratio*2^exp_N));
grid on;
