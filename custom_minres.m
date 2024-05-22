function [z, flag, relative_residual, iter, resvec] = custom_minres(A,y,eps,max_it)
    resvec = [];
    flag = 1;
    Q = y / norm(y);
    H = [];

    for i=2:max_it
        % Expand Lanczos process
        [Q, H] = lanczos(Q, H, A, i);
    
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
        c(1) = norm(y);
        Q_h_c = (Q_h')*c;
        R_h0 = R_h(1:i-1,:);
        z = R_h0 \ Q_h_c(1:i-1);
    
        % Check for convergence
        err = R_h0*z - Q_h_c(1:i-1);
        err(i) = Q_h_c(i);
        residual = norm(err(i));
        resvec = [resvec; residual];
        relative_residual = residual/norm(y);
        if relative_residual < eps
            flag = 0;
            iter = i;
            break;
        end
    end
    
end