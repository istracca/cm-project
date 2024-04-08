function [Q,R] = qr_iterative(H, Q, R, n)
    Q(n,n) = 1;
    Q_tr = Q';
    if n <= 3
        R(1:n-1,n-1) = Q_tr(1:n-1,:) * H(:,n-1);
    else
        R(n-3:n-1,n-1) = Q_tr(n-3:n-1,:) * H(:,n-1);
    end
    R(n,n-1) = H(n,n-1);
    %R = Q'*H
    c = R(n-1,n-1)/sqrt(R(n-1,n-1)^2 + R(n,n-1)^2);
    s = R(n,n-1)/sqrt(R(n-1,n-1)^2 + R(n,n-1)^2);
    Rot = eye(n);
    Rot(n-1,n-1) = c;
    Rot(n-1,n) = -s;
    Rot(n,n-1) = s;
    Rot(n,n) = c;
    R = Rot' * R;
    Q = Q * Rot;
end