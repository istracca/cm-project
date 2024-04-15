[A, y] = create_matrix_from_dimacs('netgen-3000-1-1-a-b-ns.dmx', 2312);
eps = 1e-6;

Q = y / norm(y);
H = [];
[Q, H] = lanczos(Q, H, A, 2);

c1 = H(1,1)/sqrt(H(1,1)^2 + H(2,1)^2);
s1 = H(2,1)/sqrt(H(1,1)^2 + H(2,1)^2);
Q_h = [c1 -s1; s1 c1];
R_h = Q_h' * H;


for i=3:10000
    [Q, H] = lanczos(Q, H, A, i);
    [Q_h, R_h] = qr_iterative(H, Q_h, R_h, i);
    c = zeros(i, 1);
    c(1) = norm(y);
    Q_h_c = (Q_h')*c;
    R_h0 = R_h(1:i-1,:);
    z = R_h0 \ Q_h_c(1:i-1);
    err = R_h0*z - Q_h_c(1:i-1);
    err(i) = Q_h_c(i);
    norm(err(i))
    if norm(err(i)) < eps
        fprintf("Number of iterations: %d\n", i);
        break;
    end
end