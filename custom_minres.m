function [x, flag, relative_residual, i, resvec] = custom_minres(A,y,eps,max_it)
% This function implements our custom MINRES algotirhm. 
% The choice of our input parameters is based on the MATLAB's MINRES
% function.
% Inputs:
% - A: The input matrix for which we want to solve the linear system A*x = y.
% - y: The right-hand side vector of the linear system.
% - eps: Maximum relative error to determine that convergence has occurred
% - max_it: max_it: number of iterations after which the algorithm will
%           stop, regardless of whether it has reached convergence or not.
% Returns:
% - x: solution of the linear system obtained by convergence or by
%       reaching max_it iterations
% - flag: A flag indicating whether the algorithm converged 
%       (0 if converged, 1 otherwise).
% - relative_residual: The final relative residual norm(b-A*x)/norm(b).
% - i: The number of iterations performed.
% - resvec: The vector of residuals at each iteration.


    % Initialize the residual vector
    resvec = zeros(1,max_it);
    % Initialize the flag for convergence 
    flag = 1;
    % Normalize the initial vector y
    Q = y / norm(y);
    % Initialize the Hessemberg matrix
    H = [];

    for i=2:max_it
        % Expand the Lanczos process to find the next vector and update H
        [Q, H] = lanczos(Q, H, A, i);
    
        % Calculate QR decomposition for the Hessenberg matrix
        if i==2
            % Perform initial QR decomposition for the first iteration
            c1 = H(1,1)/sqrt(H(1,1)^2 + H(2,1)^2);
            s1 = H(2,1)/sqrt(H(1,1)^2 + H(2,1)^2);
            Q_h = [c1 -s1; s1 c1];
            R_h = Q_h' * H;
        else
            % Update QR decomposition iteratively for subsequent iterations
            [Q_h, R_h] = qr_iterative(H, Q_h, R_h, i);
        end

        % Solve the linear system
        c = zeros(i, 1);
        c(1) = norm(y);
        Q_h_c = (Q_h')*c;
        R_h0 = R_h(1:i-1,:);
        z = R_h0 \ Q_h_c(1:i-1);
    
        % Check for convergence by calculating the residual
        err = R_h0*z - Q_h_c(1:i-1);
        err(i) = Q_h_c(i);
        z_ext = [z; 0];
        
        % Compute the norm of the residual
        residual = norm(err);
        % Store the residual in the vector
        resvec(i-1) = residual;
        % Calculate the relative residual
        relative_residual = residual/norm(y);

        % If the relative residual is below the tolerance, stop the iteration
        if relative_residual < eps
            % Trim the residual vector
            resvec = resvec(1:i-1);
            % Indicate convergence
            flag = 0;
            break;
        end
    end

    % Compute the final solution vector
    x = A*Q*z_ext;
end