% Implementation of the LU decomposition without a pivot.
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
    [L, U] = lu_nopivot(A);
    t2 = toc(t1);
    error = find_error(A - L*U,A);
    fprintf("[Computed in %f] Result of matrix %s : ",t2,name);
    disp(error);
end

function [L, U] = lu_nopivot(A)

n = size(A, 1); % Obtain number of rows (should equal number of columns)
L = eye(n); % L is identity matrix (for step 0, LU = (IA) = A)

for k = 1 : n
    % Remember that l_ki = a_ki / a_ii
    % So our matrix L for the one step is (L(k+1:n)) =
    % ([1,l_k+1,l_k+2...],(0,1,0,0,0,0,0),..  = A(k+1 : n,k) / A(k,k);
    leading_factor = A(k,k);

    if leading_factor == 0
        continue;
    end
    multiply_factor = A(k+1 : n,k) / leading_factor;
    L(k + 1 : n, k) = multiply_factor;

    % EXAMPLE
    % 1 -1 2
    % 4 0 1
    % -3 1 -3
    % For this matrix, the first step is dividing each row by the leading
    % number for the first row
    % for row 2 -> 4/1           1
    % for row 3 -> -3/1          4 1
    % so our L becomes          -3 0 1
    % After next step for L, we will just look at the second column. 

    % Now to step forward in calculating U, use the fact that A = LU  --> A
    % = L^-1 * A  = U for single step. The values in L give us the multiply
    % factor to make row's k'th number 0 in each row.
    for r = k + 1 : n
        deleting_row = L(r,k) * A(k,:);
        A(r, :) = A(r, :) - deleting_row; % Update each row (except the first one)
        % 1 -1 2               1 -1 2
        % 4 0 1    now becomes 4 0 3  - (4 -4 8) (multiply factor for deleting the first number * (1 -1 2 ))
        % -3 1 -3             -3 1 -3 - (3 -3 6) (multiply factor for deleting the first number * (1 -1 2 ))
    end
end
U = A; % Now A gives us the upper triangle , so it is actually U.
end


