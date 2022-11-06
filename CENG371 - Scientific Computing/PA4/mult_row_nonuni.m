function D = mult_row_nonuni(A,B,c,f)
    %{
    Input: An m × n matrix A, an n × p matrix B, a positive integer c, and probabilities {pi}
    n
    i=1.
    Output: Matrices C and R such that CR ≈ AB
    1: for t = 1 to c do
    2: Pick it ∈ {1, . . . , n} with probability Pr [it = k] = pk, in i.i.d. trials, with replacement
    3: Set C
    (t) = A(it)/
    √cpit and R(t) = B(it)/
    √cpit
    .
    4: end for
    5: Return C and R.

    %}
    [m,b] = size(A);
    [n,p] = size(B);
    C = zeros(m,n);
    R = zeros(n,p);

    p_i = ones(1,n);
    
    normalizer = 0;
    for i=1:n
        A_ri = A(:,i);
        B_ci = B(i,:);     
        p_i(1,i) = norm(A_ri)*norm(B_ci);
        normalizer = normalizer + p_i(1,i);
    end
  
    for i=1:n
        p_i(1,i) = p_i(1,i)/normalizer;
    end

    for t=1:c
        var = 1:n;
        i_t = randsample(var,1,true,p_i);
        p_it = p_i(1,i_t);

        C(:,t) = A(:,i_t)/sqrt(c*p_it);
        R(t,:) = B(i_t,:)/sqrt(c*p_it);
    end

    D = f(C,R);
       
end
