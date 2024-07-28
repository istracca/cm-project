function schur_component_trimmed = trim_schur_component(schur_component, threshold, k)
% This function trims the schur_component matrix. There are two options for trimming:
% - Strict trimming: Removes all non-diagonal elements (k = 0).
% - Soft trimming: Keeps only the k largest diagonal elements above a specified threshold ( k > 0).
% Inputs:
% - schur_component: The input Schur matrix to be trimmed.
% - threshold: The threshold value; non-diagonal elements below this value are removed.
% - k: The number of largest diagonal elements to keep.
% Returns:
% - schur_component_trimmed: The resulting trimmed Schur matrix.

    % Get the dimensions of the square matrix
    m = size(schur_component,1);

    % Create a mask for the non-diagonal elements
    off_diag = tril(ones(size(schur_component))) - eye(m);

    % Extract the non-diagonal elements
    off_diag_elements = schur_component .* off_diag;

    if k > 0
        % Find the absolute values above the threshold
        valid_elements = off_diag_elements(abs(off_diag_elements) > threshold);
    
        % Sort the valid values in descending order of their absolute values
        abs_valid_elements_sorted = sort(abs(valid_elements), 'descend');
    
        % Find the k-th largest value among the valid elements
        if length(abs_valid_elements_sorted) >= k
            kth_value = abs_valid_elements_sorted(k);
        else
            % If there are fewer than k valid elements, keep all of them
            kth_value = threshold; 
        end
        
        % Create a mask by zeroing out the non-diagonal elements above the kth_value
        off_diag_elements(abs(off_diag_elements) >= kth_value) = 0;
    end

    % Construct the final matrix by subtracting the created mask
    schur_component_trimmed = schur_component - off_diag_elements;
end
