addpath(genpath('integration'));
addpath('models');
addpath('test');
domain = [-1.5, 1.5; -1.5, 1.5];
resolution = [200,200];
initialPosition = initialize_ic_grid(resolution,domain);
n = length(initialPosition);
maxdevs = zeros(n,1);


x0 = initialPosition(10,:);
maxavg = duffing_SDE(x0, 10, 0.5)
%parfor i = 1:n
%    x0 = initialPosition(i,:);
%    maxavg = maxDevAvg(x0);
%    maxdevs(i) = maxavg;
% end
% 
% pool.delete()
% 
% 
% save('maxdevsSDE.mat','maxdevs');



