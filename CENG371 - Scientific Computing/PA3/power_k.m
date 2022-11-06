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

C = load("can_229.mat");
C = full(C.Problem.A);
currMatrix = C;

currMatrix = C; % set the current matrix for the results 

k = 229; 

power_tic = tic();
[eigVal, v] = PowerKMethod(currMatrix,k);
power_toc = toc(power_tic);
fprintf("(K = {%d})\n",k);
PrettyPrint(eigVal,v,k,power_toc,currMatrix)


%Check whether its true or not
%errorVec = CheckResult(v,eigVal,A);
%fprintf("{ %f }\n", errorVec);
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
        %fprintf('The corresponding eigenvector is: [');
        %fprintf('%g ', v(:,i));
        %fprintf(']\n\n');
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

function [eigVals,eigVecs] = PowerKMethod(A,k)
    v = ones([length(A), 1]);
    s = 1e4;    
    tolerance = 1e-14;
    eigVecs = [];
    eigVals = [];
    for j = 1:k
        for i = 1:s
            %old_norm = norm(v,1);
            v = A*v;
            %new_norm = norm(v,1);
            v = v/norm(v); % normalize predicted eigenvector
            
           % if abs(old_norm - new_norm) < tolerance
           %     break;
            %end
        end
        eigVec = v;
        eigVal = rayleigh_quotient(A,v);
        A = A - eigVal*(v*transpose(v))/(transpose(v)*v);
        eigVecs = [eigVecs,eigVec];
        eigVals = [eigVals,eigVal];
    end

end
