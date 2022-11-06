% Implementation of the LU decomposition with a pivot.
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
    [L, U, P] = lu_pivot(A);
    t2 = toc(t1);
    error = find_error(P*A - L*U,P*A);
    fprintf("[Computed in %f] Result of matrix %s : ",t2,name);
    disp(error);
end
function [L, U, P] = lu_pivot(A)

n = size(A,1);
L = eye(n); 
P = eye(n);
U = A;

for k=1:n
    %Find Pi s.t A[pi,i] != 0  --> or select the max Pi among each row
    [pv, m] = max(abs(U(k:n,k)));  % For column k,  we pick the largest row as the pivot, m is the index
    m = m + k - 1; % since we are taking from k:n , so index we found (m) is the true index m + k - 1 

    if k ~= 1 % for k=1, we dont have to change anything
        temp3=L(k,1:k-1); % do not change the L_kk ! 
        L(k,1:k-1)=L(m,1:k-1);
        L(m,1:k-1)=temp3;
    end % end of if scope

    if m~=k 

        temp1 = U(k,:); % we have to change rows in permutation matrix P
        temp2 = P(k,:); % we also have to change rows in permutation matrix P

        U(k,:) = U(m,:);
        P(k,:) = P(m,:);

        U(m,:) = temp1;
        P(m,:) = temp2;

    end

    if U(k,k) == 0
        continue;
    end
    % using the normal LU factorization algorithm
    for r=k+1:n  
        factor = U(k,k);
        L(r,k)=U(r,k)/factor;
      
        subtracted_row = L(r,k) * U(k,:); % scale by U(r,k)/factor of the row k
        U(r,:)=U(r,:)-subtracted_row; % all left is to subtract two equations
    end

    end
end

