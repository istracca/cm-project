function [x, flag, relative_residual, i, resvec] = custom_minres_preconditioned(A,y,eps,max_it,R)
% Implementation of the preconditioned version of minres. The input and 
% output parameters are the same as in the version without preconditioner, 
% except for the additional parameter R representing the matrix to be used 
% to precondition the linear system.
% The relative error for checking convergence is calculated on the original
% system, not on the preconditioned one, so that the performance of the
% preconditioned version can truly be compared with that of the
% non-preconditioned version.
% Inputs:
% - A: The input matrix for which we want to solve the linear system A*x = y.
% - y: The right-hand side vector of the linear system.
% - eps: Maximum relative error to determine that convergence has occurred
% - max_it: max_it: number of iterations after which the algorithm will
%           stop, regardless of whether it has reached convergence or not.
% - R: preconditioning matrix. The preconditioned version of the linear
%           system will be R'^(-1)AR^(-1)*Rx=R'^(-1)*y
% Returns:
% - x: solution of the original linear system obtained by convergence or by
%       reaching max_it iterations
% - flag: A flag indicating whether the algorithm converged 
%       (0 if converged, 1 otherwise).
% - relative_residual: relative error on the original system, obtained by
%       inserting the obtained x into the system
% - i: number of iterations executed by the algorithm
% - resvec: The vector of residuals at each iteration.


    % We compute only once the inverse of R
    R_inverse = (R)^(-1);

    % Compute the preconditioned version of A and y
    preconditioned_A = (R_inverse'*A*R_inverse);
    preconditioned_y = R_inverse'*y;

    % Preallocate resvec to avoid dynamic allocation
    resvec = zeros(1,max_it);
    % Initialize the flag for convergence 
    flag = 1;
    % Initialization of Q and H at iteration i=1
    Q = preconditioned_y / norm(preconditioned_y);
    H = [];

    for i=2:max_it
        % Expand the Lanczos process to find the next vector and update H
        [Q, H] = lanczos(Q, H, preconditioned_A, i);
    
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
        c(1) = norm(preconditioned_y);
        Q_h_c = (Q_h')*c;
        R_h0 = R_h(1:i-1,:);
        z = R_h0 \ Q_h_c(1:i-1);
    
        % Obtain x from z
        z_ext = [z; 0];
        preconditioned_x = Q*z_ext;
        x = R_inverse*preconditioned_x;
        
        % Compute the error on the original system and update of resvec
        err = A*x - y;
        % Compute the norm of the residual
        residual = norm(err);
        % Store the residual in the vector
        resvec(i-1) = residual;
        % Calculate the relative residual
        relative_residual = residual/norm(y);

        % Check for convergence
        if relative_residual < eps
            % Trim the residual vector
            resvec = resvec(1:i-1);
            % Indicate convergence
            flag = 0;
            break;
        end
    end
    
end