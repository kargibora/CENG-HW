
function C = mult_naive(A,B)
    % Algorithm 2
    %{
    for i = 1 to m do
    2: for j = 1 to p do
    3: (AB)ij = 0
    4: for k = 1 to n do
    5: (AB)ik += AijBjk    
    6: end for
    7: end for
    8: end for
    9: Return AB
    %}
    [m,n] = size(A);
    [n,p] = size(B);
    C = zeros(m,p);
    for i=1:m
        for j=1:n
            for k=1:p
                C(i,k) = C(i,k) + A(i,j)*B(j,k);
            end
        end
    end

end