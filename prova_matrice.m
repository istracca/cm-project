[A, y] = create_matrix_from_dimacs('graph_instances/net12_8_1.dmx', 1);
prova = 1
% Error threshold
eps = 1e-6;

% Misura tempo per la prima sezione
tic;

[custom_x, custom_flag, custom_relres, custom_iter, custom_resvec] = custom_minres(A, y, 1e-6, 10000);
custom_relative_resvec = custom_resvec/norm(y);

elapsedTime1 = toc;
disp(['Tempo per la creazione della matrice: ', num2str(elapsedTime1), ' secondi']);


% Misura tempo per la seconda sezione
tic;
[x, flag, relres, iter, resvec] = minres(A, y, 1e-6, 10000);
relative_resvec = resvec/norm(y);
%relative_resvec = relative_resvec(2:end);

elapsedTime2 = toc;
disp(['Tempo per il calcolo della matrice al quadrato: ', num2str(elapsedTime2), ' secondi']);

% Grafico dell'evoluzione dei residui relativi
figure;
semilogy(custom_relative_resvec, 'b', 'LineWidth', 1.5);
hold on;
semilogy(relative_resvec, 'r', 'LineWidth', 1.5);
hold off;

xlabel('Iterazioni');
ylabel('Residuo relativo');
legend('Custom MINRES', 'Matlab MINRES');
title('Evoluzione dei residui relativi');
grid on;
