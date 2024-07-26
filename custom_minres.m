% Implementazione del nostro minres. Nella scelta dei parametri di input e
% di output ci siamo ispirati al minres di matlab


function [x, flag, relative_residual, i, resvec] = custom_minres(A,y,eps,max_it)
    resvec = zeros(1,max_it);
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
        err = R_h0*z - Q_h_c(1:i-1);        % should be zero
        err(i) = Q_h_c(i);
        z_ext = [z; 0];
           
        residual = norm(err);
        resvec(i-1) = residual;
        relative_residual = residual/norm(y);
        if relative_residual < eps
            resvec = resvec(1:i-1);
            flag = 0;
            break;
        end
    end
    x = A*Q*z_ext;
end