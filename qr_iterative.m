% Function that implements Givens rotations to obtain QR decomposizion of H 
% Input parameters:
% - H: The new tridiagonal matrix (dimensions nx(n-1)) that is being
%       factorized, given the factorization of the old one
% - Q: The orthogonal matrix (dimensions (n-1)x(n-1)) that is part of the
%       QR decomposition of the old H (that was ((n-1)x(n-2))
% - R: The upper triangular matrix (dimensions ((n-1)x(n-2)) that is part
%       of the QR decomposition of the old H (that was ((n-1)x(n-2))
% - n: the current step or iteration number in the process that updates Q 
%      and R.
% Returns:
% - Q: The orthogonal matrix (dimensions nxn) that is part of the QR
%       decomposition of the new H
% - R: The upper triangular matrix (dimensions (nx(n-1)) that is part of
%       the QR decopmosition of the new H 
function [Q,R] = qr_iterative(H, Q, R, n)
    % Additional row and column for Q
    Q(n,n) = 1;
    
    Q_tr = Q';

    % R update 
    if n <= 3
        % The whole column is updated
        R(1:n-1,n-1) = Q_tr(1:n-1,:) * H(:,n-1);
    else
        % Update only the last three elements of the column
        R(n-3:n-1,n-1) = Q_tr(n-3:n-1,:) * H(:,n-1);
    end

    R(n,n-1) = H(n,n-1);

    % Cosine and sine update
    c = R(n-1,n-1)/sqrt(R(n-1,n-1)^2 + R(n,n-1)^2);
    s = R(n,n-1)/sqrt(R(n-1,n-1)^2 + R(n,n-1)^2);
    
    % Rotation matrix construction
    Rot = eye(n);
    Rot(n-1,n-1) = c;
    Rot(n-1,n) = -s;
    Rot(n,n-1) = s;
    Rot(n,n) = c;

    % Givens rotation
    R = Rot' * R;
    Q = Q * Rot;
end