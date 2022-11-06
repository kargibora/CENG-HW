% Implementation of inverse power method
A = [2 -1  0  0  0;
    -1  2 -1  0  0;
     0 -1  2 -1  0;
     0  0 -1  2 -1;
     0  0  0 -1  2];

B = [0.2  0.3 -0.5;
     0.6 -0.8  0.2;
    -1.0, 0.1  0.9];

alpha = 0;
currMatrix = A;

ipm_tic = tic();
[eigVal, v] = InversePowerMethod(currMatrix, alpha);
ipm_time = toc(ipm_tic);

PrettyPrint(eigVal,v,1,ipm_time,currMatrix,alpha)
%Check whether its true or not
%errorVec = CheckResult(v,eigVal,A);
%fprintf("{ %f }\n", errorVec);
function error = CheckResult(eigVec,eigVal,A)
    result = A * eigVec;
    eigResult = eigVal*eigVec;
    error = eigResult - result; % this should be close to 0 if everything is true
end

function PrettyPrint(eigVal,v,k,t,A,alpha)
% function for printing the results clearly 
fprintf("Total runtime : %f \n",t);
    for i=1:k
        fprintf("The eigenvalue closes to %f is %f \n",alpha,eigVal(i));
        fprintf('The corresponding eigenvector is: [');
        fprintf('%g ', v(:,i));
        fprintf(']\n\n');
        %error = CheckResult(v(:,i),eigVal(i),A);
        %fprintf(" Total error = %f",sum(error));
    end
end



function [eigVal] = RayleighQuotient(A,v)
    % if v is a good approximation, rayleigh quotient method
    % gives us a good approximation for the eigenvalue
    vt = transpose(v);
    vtv = vt*v;
    Av = A*v;
    vtav = vt*Av;
    eigVal = vtav/vtv;
end

function [eigVal,eigVec] = InversePowerMethod(A,alpha)
    alphaI = alpha*eye(length(A));
    v = ones([length(A),1]);
    v = v / norm(v);
    A_shifted = A - alphaI;

    k = 1e7;
    tolerance = 1e-5;
    A_shifted_inverse = inv(A_shifted);
    for i = 1:k
        prev_norm = norm(v);
        v = A_shifted_inverse*v;

        if abs(prev_norm - norm(v)) < tolerance % check for a speedup
            break;
        end

        v = v/norm(v); % normalize predicted eigenvector
    end
    eigVec = v;
    eigVal = RayleighQuotient(A,v); % we found lambda close to alpha
end

