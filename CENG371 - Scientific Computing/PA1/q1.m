N = 1000;

Lis = [];
Res = [];
Zeros = [];
for k=1:N % from 1 to N
	sum = k * ((k + 1)/k - 1) - 1 ;
    Lis = [Lis ; k];
    Res = [Res ; sum/eps];
    if sum/eps == 0
        Zeros = [Zeros ; k];
    end
      
end
plot(Lis,Res);

fprintf('Zero values: [%s]\n', join(string(Zeros), ','));


    