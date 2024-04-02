function [Q,R] = qr_iterative(H, Q, R, n)
    Q(n,n) = 1;
    R = Q'*H;
    c = H(n-1,n-1)/sqrt(H(n-1,n-1)^2 + H(n,n-1)^2);
    s = H(n,n-1)/sqrt(H(n-1,n-1)^2 + H(n,n-1)^2);
    Rot = eye(n);
    Rot(n-1,n-1) = c;
    Rot(n-1,n) = -s;
    Rot(n,n-1) = s;
    Rot(n,n) = c;
    R = Rot' * R;
    Q = Q * Rot;
end