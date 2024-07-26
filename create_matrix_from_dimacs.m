% Ha senso dividere e fare un file per ogni funzione

function [D, E, b] = create_matrix_from_dimacs(file, seed)
    % Ottieni il grafico dal file
    graph = create_graph_from_dimacs(file);
    % Ottieni la matrice di incidenza (completa) dal grafico
    E = incidence_matrix(graph);
    
    m = size(E, 2);
    n = graph.num_nodes;
    E = E(1:n-1, :);
    
    rng(seed);

    % decommenta questa parte per generare D di tutti uni
    % D_values = ones(m,1);

    % decommenta questa parte per elementi di D in [0,1]
    % D_values = rand(m,1);
    % D_values = D_values*19 + 1;
    

    % decommenta questa parte per elementi di D in cluster [0.95,1.05] e
    % [19.95,20.05
    half_m = floor(m / 2);
    D_cluster1 = 0.95 + (0.1) * rand(half_m, 1);
    D_cluster2 = 19.95 + (0.1) * rand(m - half_m, 1);
    D_values = [D_cluster1; D_cluster2];
    D_values = D_values(randperm(m));

    % Costruisci D partendo da D_values
    D = spdiags(D_values, 0, m, m);
    % Ottieni b come concatenazione di costi degli archi e flussi dei nodi
    % del grafo
    b = costs_flows(graph);
end

function graph = create_graph_from_dimacs(file)
    graph = struct();
    graph.edges = [];
    graph.sources_sinks = [];

    fid = fopen(file, 'r');
    if fid == -1
        error('Error in file opening\n');
    end
    count = 0;
    line = fgetl(fid);

    while ~(line == -1)
        tokens = strsplit(line);
        if tokens{1} == 'c'
            count = count + 1;
            if count == 5
                line = strrep(line, ' ', '');
                line = split(line, ':');
                num_nodes = line(2);
                num_nodes = str2num(num_nodes{1});
            end
        elseif tokens{1} == 'a'
            edge = str2double(tokens(2:end));
            graph.edges(end+1, :) = edge;
        elseif tokens{1} == 'n'
            node = str2double(tokens(2:end));
            graph.sources_sinks(end+1, :) = node;
        end
        line = fgetl(fid);
    end
    graph.num_nodes = num_nodes;
    fclose(fid);
end

function incidence_matrix = incidence_matrix(graph)
    edges = graph.edges;
    num_nodes = graph.num_nodes;
    num_edges = size(edges, 1);
    incidence_matrix = zeros(num_nodes, num_edges);

    for i = 1:num_edges
        source = edges(i,1);
        dest = edges(i,2);
        incidence_matrix(source, i) = 1;
        incidence_matrix(dest, i) = -1;
    end
end


function costs_flows = costs_flows(graph)
    flows = zeros(graph.num_nodes-1,1);
    for i=1:length(graph.sources_sinks)-1
        node_index = graph.sources_sinks(i,1);
        flows(node_index) = graph.sources_sinks(i,2);
    end
    costs = graph.edges(:,4);
    costs_flows = [-costs;flows];
end