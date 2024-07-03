function [x, flag, relative_residual, i, resvec] = custom_minres_preconditioned(A,y,eps,max_it,P)
    preconditioned_A = P*A*P';
    preconditioned_y = P*y;
    resvec = [];
    Q = preconditioned_y / norm(preconditioned_y);
    H = [];
    flag = 1;

    for i=2:max_it
        % Expand Lanczos process
        [Q, H] = lanczos(Q, H, preconditioned_A, i);
    
        % Calculate QR decomposition
        if i==2
            c1 = H(1,1)/sqrt(H(1,1)^2 + H(2,1)^2);
            s1 = H(2,1)/sqrt(H(1,1)^2 + H(2,1)^2);
            Q_h = [c1 -s1; s1 c1];
            R_h = Q_h' * H;
        else
            [Q_h, R_h] = qr_iterative(H, Q_h, R_h, i);
        end
        % Solving the linear system
        c = zeros(i, 1);
        c(1) = norm(preconditioned_y);
        Q_h_c = (Q_h')*c;
        R_h0 = R_h(1:i-1,:);
        z = R_h0 \ Q_h_c(1:i-1);
    
        % Check for convergence
        z_ext = [z; 0];
        preconditioned_x = Q*z_ext;
        x = P'*preconditioned_x;
        err = A*x - y;
        residual = norm(err);
        resvec = [resvec; residual];
        relative_residual = residual/norm(y);
        if relative_residual < eps
            flag = 0;
            break;
        end
    end
    
end