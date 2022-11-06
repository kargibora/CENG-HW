function D = mult_proj_Gauss_orth(A,B,c,f)
    [m,n] = size(A);
    [n,p] = size(B);
    P = normrnd(0,1/sqrt(n),n,c);

    % Gram schmid P = [n,c]
    [Q,R_] = qr(P);
    C = f(A,Q);
    R = f(transpose(Q),B);
    D = f(C,R);
end