addpath(genpath('integration'));
addpath('models');
addpath('test');
domain = [-1.5, 1.5; -1.5, 1.5];
resolution = [100,100];
initialPosition = initialize_ic_grid(resolution,domain);
n = length(initialPosition);
maxdevs = zeros(n,1);
addpath('models');
addpath('integration');
addpath('integration/grid');
euler = parcluster('local');
pool = parpool(euler,8);

parfor i = 1:n
    x0 = initialPosition(i,:);
    maxdevs(i) = maxDev(@(t,x,e) d_hopf(t,x,e),x0, 0:0.001*2*pi:3*2*pi);
end

pool.delete()


save('maxdevsHopf.mat','maxdevs');
