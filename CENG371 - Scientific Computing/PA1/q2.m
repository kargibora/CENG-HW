N = 1e6;
fsum = 0;
number_arange = [1:10^6]; 
Nums = 1 + (1e6 + 1 - number_arange)*1e-8;
single_Nums = single(Nums);
theoritical_value = 1005000.005;

time_naive = tic(); % Setup naive sum
result_naive = naive_sum(Nums,N);
time_naive = toc(time_naive);

time_pairwise = tic(); % Setup pairwise sum
result_pairwise= pairwise_sum(Nums,N);
time_pairwise = toc(time_pairwise);

time_compansated = tic(); % Setup compansated sum
result_compansated = kahan_compansated_sum(Nums,N);
time_compansated = toc(time_compansated);

%fprintf("Naive Sum : %.16f   Time : %.16f  Error : %.16f  Relative Error : %.16f \n",result_naive,time_naive, result_naive - theoritical_value, (result_naive - theoritical_value)/theoritical_value);
%fprintf("Pairwise Sum:%.16f   Time : %.16f Error : %.16f  Relative Error : %.16f \n",result_pairwise,time_pairwise,result_pairwise - theoritical_value, (result_pairwise - theoritical_value)/theoritical_value);
%fprintf("Compansated Sum : %.16f    Time : %.16f Error : %.16f  Relative Error : %.16f \n",result_compansated,time_compansated,result_compansated - result_compansated, (result_compansated - theoritical_value)/theoritical_value);
%fprintf("Pairwise sum algorithm for base case N = 10000\n")
pretty_print(result_naive,time_naive,theoritical_value,"NAIVE SUM ALGORTIHM");
pretty_print(result_pairwise,time_pairwise,theoritical_value,"PAIRWISE SUM ALGORTIHM");
pretty_print(result_compansated,time_compansated,theoritical_value,"COMPANSATED SUM ALGORTIHM");


time_naive = tic(); % Setup naive sum
result_naive = naive_sum(single_Nums,N);
time_naive = toc(time_naive);

time_pairwise = tic(); % Setup pairwise sum
result_pairwise= pairwise_sum(single_Nums,N);
time_pairwise = toc(time_pairwise);

time_compansated = tic(); % Setup compansated sum
result_compansated = kahan_compansated_sum(single_Nums,N);
time_compansated = toc(time_compansated);

%fprintf("Naive Sum : %.16f   Time : %.16f  Error : %.16f  Relative Error : %.16f \n",result_naive,time_naive, result_naive - theoritical_value, (result_naive - theoritical_value)/theoritical_value);
%fprintf("Pairwise Sum:%.16f   Time : %.16f Error : %.16f  Relative Error : %.16f \n",result_pairwise,time_pairwise,result_pairwise - theoritical_value, (result_pairwise - theoritical_value)/theoritical_value);
%fprintf("Compansated Sum : %.16f    Time : %.16f Error : %.16f  Relative Error : %.16f \n",result_compansated,time_compansated,result_compansated - result_compansated, (result_compansated - theoritical_value)/theoritical_value);

pretty_print(result_naive,time_naive,theoritical_value,"NAIVE SUM ALGORTIHM (SINGLE)");
pretty_print(result_pairwise,time_pairwise,theoritical_value,"PAIRWISE SUM ALGORTIHM (SINGLE)");
pretty_print(result_compansated,time_compansated,theoritical_value,"COMPANSATED SUM ALGORTIHM (SINGLE)");


function pretty_print(result,time,theo,name)
fprintf("_________________________________________\n")
fprintf("ALGORTIHM : %s \n",name);
fprintf("Calculated sum : %.16f\n",result);
fprintf("Error : %.16f \n",theo - result);
fprintf("Relative Error : %.16f \n", (theo - result)/theo);
fprintf("Total execution time : %.16f\n",time)
fprintf("_________________________________________\n")
end
function s = naive_sum(Numbers,N)

    s = 0.0;
    for i=1:N
        s = s + Numbers(i);
    end
end

       
function s = kahan_compansated_sum(Numbers,N)
    s = 0.0; % Accumulator
    c = 0.0; % Compensation

    for i=1:N
        y = Numbers(i) - c;
        t = s + y;
        c = (t - s) - y;
        s = t;

    end



end




function s = pairwise_sum(Numbers,n)
base_case = 2;
   if n == 0 

       s = 0;

elseif n <= base_case
     s = 0;
     for i=1:n
         s = s + Numbers(i);
     end
 
   else
     m = floor(n/2);
     s = pairwise_sum(Numbers(1:m),m) + pairwise_sum(Numbers(m+1:n),n-m);
   end
end