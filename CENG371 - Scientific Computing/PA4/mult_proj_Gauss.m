function D = mult_proj_Gauss(A,B,c,f)
    [m,n] = size(A);
    [n,p] = size(B);
    P = normrnd(0,1/sqrt(n),n,c);
    C = f(A,P);
    R = f(transpose(P),B);
    D = f(C,R);
end