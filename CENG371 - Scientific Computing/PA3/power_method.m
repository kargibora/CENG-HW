% Implementation of Power Method
% Bora Kargı
A = [2 -1  0  0  0;
    -1  2 -1  0  0;
     0 -1  2 -1  0;
     0  0 -1  2 -1;
     0  0  0 -1  2];

B = [0.2  0.3 -0.5;
     0.6 -0.8  0.2;
    -1.0, 0.1  0.9];

currMatrix = B;
x = ones([length(currMatrix), 1]);

power_tic = tic();
[eigVal, v] = PowerMethod(currMatrix, x);
power_time = toc(power_tic);

PrettyPrint(eigVal,v,1,power_time,currMatrix)

function error = CheckResult(eigVec,eigVal,A)
    result = A * eigVec;
    eigResult = eigVal*eigVec;
    error = eigResult - result; % this should be close to 0 if everything is true
end

function PrettyPrint(eigVal,v,k,t,A)
% function for printing the results clearly 
fprintf("Total runtime : %f \n",t);
    for i=1:k
        fprintf("The %d'th largest eigenvalue is %f \n",i,eigVal(i));
        fprintf('The corresponding eigenvector is: [');
        fprintf('%g ', v(:,i));
        fprintf(']\n\n');
        %error = CheckResult(v(:,i),eigVal(i),A);
        %fprintf(" Total error = %f",sum(error));
    end
end

function [eigVal] = rayleigh_quotient(A,v)
    % if v is a good approximation, rayleigh quotient method
    % gives us a good approximation for the eigenvalue
    vt = transpose(v);
    vtv = vt*v;
    Av = A*v;
    vtav = vt*Av;
    eigVal = vtav/vtv;
end

function [eigVal,eigVec] = PowerMethod(A,v)
    k = 1e7;
    tolerance = 1e-7;
    for i = 1:k
        prev_norm = norm(v);
        v = A*v;

        if abs(prev_norm - norm(v)) < tolerance % check for a speedup
            break;
        end

        v = v/norm(v); % normalize predicted eigenvector
    end
    eigVec = v;
    eigVal = rayleigh_quotient(A,v);
end
