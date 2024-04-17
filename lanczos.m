% Function that implements Lanczos algorithm
% Input parameters:
% - Q: Orthogonal matrix whose columns form the basis of the Krylov 
%      subspace of dimension (n-1)
% - H: Tridiagonal matrix of dimensions (n-1)x(n-2)
% - A: The original matrix.
% - n: The current iteration in the Lanczos process
% Returns:
% - Q: Updated orthogonal matrix (nxn) with columns representing the basis of the
%      expanded Krylov subspace (dimension n) including the latest vector.
% - H: Updated tridiagonal matrix reflecting the orthogonal projection of
%      matrix A onto the new basis represented by Q (dimensions nx(n-1))
function [Q, H] = lanczos(Q, H, A, n)
    % First iteration is different from the following
    if n==2
        % Projection
        w = A*Q(1:end,1);
        % Orthogonalization
        H(1,1) = Q(:, 1)' * w;
        % Normalization
        w = w-Q(1:end, 1) * H(1,1);
        H(2,1) = norm(w);
        % Q update
        Q(1:end, 2) = w / H(2,1);

    % General iteration
    else
        % Projection
        w = A*Q(1:end,n-1);
        % Copy previous off-diagonal element
        H(n-2, n-1) = H(n-1, n-2);
        % Orthogonalization
        w = w - Q(1:end, n-2) * H(n-2, n-1);
        H(n-1, n-1) = Q(1:end, n-1)' * w;
        w = w-Q(1:end, n-1) * H(n-1, n-1);
        % Normalization
        H(n, n-1) = norm(w);
        Q(1:end, n) = w / norm(w);
    end
end