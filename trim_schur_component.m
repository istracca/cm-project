function schur_component_trimmed = trim_schur_component(schur_component, threshold, k)
    % Dimensioni della matrice
    m = size(schur_component,1);

    % Creare una maschera per gli elementi non diagonali inferiori
    off_diag = tril(ones(size(schur_component))) - eye(m);

    % Estrazione degli elementi non diagonali
    off_diag_elements = schur_component .* off_diag;

    if k > 0

        % Trova i valori sopra la soglia
        valid_elements = off_diag_elements(abs(off_diag_elements) > threshold);
    
        % Ordina i valori validi in ordine decrescente
        abs_valid_elements_sorted = sort(abs(valid_elements), 'descend');
    
        % Trova il k-esimo valore più grande sopra la soglia
    
        if length(abs_valid_elements_sorted) >= k
            kth_value = abs_valid_elements_sorted(k);
        else
            kth_value = threshold; % Se ci sono meno di k elementi, mantieni tutti quelli validi
        end
    
        fprintf('%.4g', kth_value);
    
        % Mantieni solo i k elementi più grandi in modulo sopra la soglia
        off_diag_elements(abs(off_diag_elements) >= kth_value) = 0;
    end

    % Costruisci la matrice finale
    schur_component_trimmed = schur_component - off_diag_elements;
end
