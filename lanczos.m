function [Q, H] = lanczos(Q, H, A, n, y)
    if n==1
        Q = y / norm(y);
        H = eye(1);

    elseif n==2
        w = A*Q(1:end,1);
        H(1,1) = Q(:, 1)' * w;
        w = w-Q(1:end, 1) * H(1,1);
        H(2,1) = norm(w);
        Q(1:end, 2) = w / H(2,1);
    
    else
        w = A*Q(1:end,n-1);
        H(n-2, n-1) = H(n-1, n-2);
        w = w - Q(1:end, n-2) * H(n-2, n-1);
        H(n-1, n-1) = Q(1:end, n-1)' * w;
        w = w-Q(1:end, n-1) * H(n-1, n-1);
        H(n, n-1) = norm(w);
        Q(1:end, n) = w / norm(w);

    end
end



