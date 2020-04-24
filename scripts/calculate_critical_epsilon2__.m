addpath(genpath('integration'));
addpath('models');
addpath('test');
lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
resolution = [100,100];
domain = [-1.5, 1.5; -1.5, 1.5];

initialPosition = initialize_ic_grid(resolution,domain);
euler = parcluster('local');
pool = parpool(euler,32);

epsilons = logspace(1, 3, 10000);
criteps2 = zeros(length(initialPosition), 1);
T = 3*2*pi;
deltaT = T/15;
timeInterval = [0, T];
stepSize = deltaT;
parfor i = 1:length(initialPosition)
   x0 = initialPosition(i,:);
   lDerivative = @(t,y) d_phi(t,y,0, false);
   [t, uncertAtTime] = computeAlongTrajectory(lDerivative, x0, timeInterval, stepSize );
   lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
   criteps2(i) = getCriticalEpsilon(lDerivative,x0,uncertAtTime, t, epsilons);
end

pool.delete();
save('criticalEps__FFINE.mat', 'criteps2');
