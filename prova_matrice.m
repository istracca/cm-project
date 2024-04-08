[E, D, y] = utility_read_matrix('netgen-3000-1-1-a-b-ns.dmx', 2312, 0);
D = diag(D);
A = [D E'; E zeros(154)];

Q = y / norm(y);
H = [];

[Q, H] = lanczos(Q, H, A, 2);

c1 = H(1,1)/sqrt(H(1,1)^2 + H(2,1)^2);
s1 = H(2,1)/sqrt(H(1,1)^2 + H(2,1)^2);
Q_h = [c1 -s1; s1 c1];

R_h = Q_h' * H

for i=3:10
    [Q, H] = arnoldi(Q, H, A, i);
    [Q_h, R_h] = qr_iterative(H, Q_h, R_h, i);
        

end