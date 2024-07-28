function [D, E, b] = create_matrix_from_dimacs(file, flag, seed)
% This function reads a graph from a DIMACS file and returns a sparse
% diagonal matrix D, the incidence matrix without the last row E and a 
% vector b. The matrix [D E; E' 0] and the vector b are typically used in
% flow network optimization problems.
% The function first constructs a graph structure, then derives an
% incidence matrix from the graph. It uses the incidence matrix and
% randomly generated diagonal matrix D to form the matrix A. The vector b
% is generated based on the costs and flows defined in the graph.
% Input parameters:
% - file: DIMACS file 
% - flag: a flag to determine the policy of construction of D. It can be:
%       - 0: D values all equal to one
%       - 1: D values randomly uniformly distributed in [0,1]
%       - 2: D values randomly uniformly distributed in [1,20]
%       - 3: D values randomly uniformly distributed in two clusters:
%       [0.95,1.05] and [19.95,20.05]
% - seed: random seed 
% Returns:
% - D: sparse diagonal matrix. Its construction depends on the input flag
% - b: vector containing combined costs and flows.

    % Constuct the graph from DIMACS file
    graph = create_graph_from_dimacs(file);

    % Create the incidence matrix from the graph
    E = incidence_matrix(graph);
    
    % Determine the number of edges and nodes in the graph
    m = size(E, 2);
    n = graph.num_nodes;

    % Use only the first n-1 rows of E
    E = E(1:n-1, :);
    
    % Fix the random seed to allow replicability (it is needed in the 
    % generation of D_values).
    rng(seed);

    % Generate D_values according to input flag
    switch flag
        case 0
            D_values = ones(m,1);
        case 1
            D_values = rand(m,1);
        case 2
            D_values = rand(m,1);
            D_values = D_values*19 + 1;
        case 3
            half_m = floor(m / 2);
            D_cluster1 = 0.95 + (0.1) * rand(half_m, 1);
            D_cluster2 = 19.95 + (0.1) * rand(m - half_m, 1);
            D_values = [D_cluster1; D_cluster2];
            D_values = D_values(randperm(m));
    end

    % Create a sparse diagonal D from D_values
    D = spdiags(D_values, 0, m, m);
    % Get b as concatenation of graph edge costs and graph node flows
    b = costs_flows(graph);
end

function graph = create_graph_from_dimacs(file)
% This function reads a graph from a DIMACS formatted file and constructs
% a graph data structure with edges, sources, sinks, and the total number of nodes.
% Inputs:
% - file: DIMACS file.
% Outputs:
% - graph: structure containing the following fields:
%   - edges: matrix where each row represents an edge with
%       start node, end node, capacity, and optionally a cost.
%   - sources_sinks: matrix where each row represents a source or sink
%       with node number and associated flow value.
%   - num_nodes: integer specifying the total number of nodes in the graph.

    % Initialize an empty graph structure
    graph = struct();
    graph.edges = [];
    graph.sources_sinks = [];

    % Open the DIMACS file
    fid = fopen(file, 'r');
    if fid == -1
        error('Error in file opening\n');
    end

    count = 0;
    line = fgetl(fid);

    % Parse the file line by line
    while ~(line == -1)
        tokens = strsplit(line);
        if tokens{1} == 'c'
            % Comment counter
            count = count + 1;
            % At the fifth comment line, get the number of nodes
            if count == 5
                line = strrep(line, ' ', '');
                line = split(line, ':');
                num_nodes = line(2);
                num_nodes = str2num(num_nodes{1});
            end
        % Get the information of each arch 
        elseif tokens{1} == 'a'
            edge = str2double(tokens(2:end));
            graph.edges(end+1, :) = edge;
        % Get the information of each node (source/sink)
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
% This function generates an incidence matrix for a given graph, which
% is represented as a structure containing the nodes and edges of the graph.
% Input:
% - graph: a graph structure
% Returns:
% - incidence matrix: matrix that represents the incidence relationships between 
%   nodes and edges

    % Extract edges and number of nodes from the graph structure
    edges = graph.edges;
    num_nodes = graph.num_nodes;

    % Determine the number of edges
    num_edges = size(edges, 1);

    % Initialize the incidence matrix with zeros
    incidence_matrix = zeros(num_nodes, num_edges);

    % Iterate over each edge to populate the incidence matrix
    for i = 1:num_edges
        % Extract the source and destination nodes for the current edge
        source = edges(i,1);
        dest = edges(i,2);

        % Set the matrix entries: +1 for source, -1 for destination
        incidence_matrix(source, i) = 1;
        incidence_matrix(dest, i) = -1;
    end
end


function costs_flows = costs_flows(graph)
% This function extracts and combines edge costs and node flows from a
% graph structure into a single vector.
% Input:
% - graph: a graph structure
% Returns:
% - costs_flows: a vector combining edge costs and node flows, where the
%   edge costs are listed first followed by the node flows.

    % Initialize a vector of zeros for node flows
    flows = zeros(graph.num_nodes-1,1);

    % Assign flows based on the sources_sinks information, excluding the last entry
    for i=1:length(graph.sources_sinks)-1
        node_index = graph.sources_sinks(i,1);
        flows(node_index) = graph.sources_sinks(i,2);
    end

    % Extract costs from the fourth column of the edges array
    costs = graph.edges(:,4);

    % Combine costs and flows into a single vector    
    costs_flows = [costs;flows];
end