% Implementazione della versione precondizionata del minres. I parametri di
% input e di output sono gli stessi della versione senza preconditioner,
% fatta eccezione per il parametro R che rappresenta la matrice da usare
% per precondizionare il sistema lineare.
% l'errore relativo per controllare la convergenza viene calcolato sul
% sistema originale, e non su quello precondizionato, in modo da poter
% confrontare effettivamente le prestazioni della versione precondizionata
% cone quelle della versione non precondizionata



function [x, flag, relative_residual, i, resvec] = custom_minres_preconditioned(A,y,eps,max_it,R)
    R_inverse = (R)^(-1);
    preconditioned_A = (R_inverse'*A*R_inverse);
    preconditioned_y = R_inverse'*y;
    resvec = zeros(1,max_it);
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
        x = R_inverse*preconditioned_x;
        err = A*x - y;
        residual = norm(err);
        resvec(i-1) = residual;
        relative_residual = residual/norm(y);
        if relative_residual < eps
            resvec = resvec(1:i-1);
            flag = 0;
            break;
        end
    end
    
end