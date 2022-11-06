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


k = 229; 

sub_tic = tic();
[eigVal, v] = SubspaceIterations(currMatrix,k);
sub_time = toc(sub_tic);
fprintf("(K = {%d})\n",k);
PrettyPrint(eigVal,v,k,sub_time,currMatrix)

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
        %fprintf(" Total error = %f \n",sum(error)); can be used to check
        % errors
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

function [eigVals,eigVecs] = SubspaceIterations(A,k)
    s = 1e4;    
    tolerance = 1e-14;
    eigVecs = [];
    eigVals = [];
    X = randn([length(A),k]);

    for j = 1:s
        %old_norm = norm(X,1);
        Z = A*X;
        [X,R] = qr(Z);
        %new_norm = norm(X,1);
        %if abs(old_norm - new_norm) < tolerance
        %    break;
        %end
    end

    eigVecs = X;
    
    for i=1:k
        eigVals = [eigVals, rayleigh_quotient(A,X(:,i))];
    end
end
