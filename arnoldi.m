function [Q, H] = arnoldi(Q, H, A, n)

    if n==2
        w = A*Q(1:end,1);
        H(1,1) = Q(:, 1)' * w;
        w = w-Q(1:end, 1) * H(1,1);
        H(2,1) = norm(w);
        Q(1:end, 2) = w / H(2,1);

    else
        w = A*Q(1:end,n-1);
        for i=1:n-1
            H(i,n-1) = Q(1:end,i)' * w;
            w = w-Q(1:end,i) * H(i,n-1);
        end
        H(n,n-1) = norm(w);
        Q(1:end,n) = w/norm(w);
    end
end



