% Implementation of the LU decomposition with a complete pivot.
% Bora KARGI
hilb_n = 1000;
A = [1,2,3,4;2,4,6,8;3,6,9,12;4,8,12,16];
B = hilb(hilb_n);
C = B;
C(1:hilb_n+1:end) = B(1:hilb_n+1:end) * 1e-16;

%%%%%%%%%%%%
%Additional%
%D = randn(1000,1000);
%calculate_result(D,"[D]");

calculate_result(A,"[A]");
calculate_result(B,"[B]");
calculate_result(C,"[C]");


function relative_error = find_error(A,B)
    norm_A = norm(A,2);
    norm_B = norm(B,2);
    relative_error = norm_A / norm_B;
end

function error = calculate_result(A,name)
    t1 = tic();
    [L, U, P, Q] = lu_pivotr(A);
    t2 = toc(t1);
    error = find_error(P*A*Q - L*U,P*A*Q);
    fprintf("[Computed in %f] Result of matrix %s : ",t2,name);
    disp(error);
end


function [L, U, P, Q] = lu_pivotr(A)
% P * A* Q = L * U
% So we permute + rotate at the same time.
n = size(A,1);

P = eye(n);
Q = eye(n);
L = eye(n);

for k = 1 : n-1
    % We need maximum leading pivot and the maximum leading column
    % A[i, 1:n] = A[pi,1:n]
    % A[1:n,i] = A[1:n,qi]
    [max_row, ris] = max(abs( A(k:n,k:n)) ); % row indexes of the maximum value in columns
    [max_col, ci] = max(max_row); % column index of the maximum value in array
%    L=zeros(n);


    ri = ris(ci); % ci is the index of maximum value among all maximum values in columns,
    % so ri is the corresponding index for the maximum value's row

    % refactor the indexes so it matches the global position in the array
    % (add k - 1)
    ri = ri + k - 1;
    ci = ci + k - 1;

    % swap row and column
    temp_row = A(k,:);
    A(k,:) = A(ri,:);
    A(ri,:) = temp_row;

    temp_col = A(:,k);
    A(:,k) = A(:,ci);
    A(:,ci) = temp_col;

    % swap row and column of P and Q accordingly
    temp_row = P(k,:);
    P(k,:) = P(ri,:);
    P(ri,:) = temp_row;

    temp_col = Q(:,k);
    Q(:,k) = Q(:,ci);
    Q(:,ci) = temp_col;

    if k ~= 1 % for k=1, we dont have to change anything
        temp3=L(k,1:k-1); % do not change the L_kk ! 
        L(k,1:k-1)=L(ri,1:k-1);
        L(ri,1:k-1)=temp3;

        temp3=L(1:k-1,k); % do not change the L_kk ! 
        L(1:k-1,k)=L(1:k-1,ci);
        L(1:k-1,ci)=temp3;

    end 

    if A(k,k) == 0
        continue;
    end
    for r=k+1:n  
        factor = A(k,k);
        L(r,k)=A(r,k)/factor;

        subtracted_row = L(r,k) * A(k,:); % scale by U(r,k)/factor of the row k
        A(r,:)=A(r,:)-subtracted_row; % all left is to subtract two equations
    end

end

U = A;
end