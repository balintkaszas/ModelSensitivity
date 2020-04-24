addpath(genpath('integration'));
addpath('models');
addpath('test');
lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
resolution = [1000,1000];
domain = [-1.5, 1.5; -1.5, 1.5];

initialPosition = initialize_ic_grid(resolution,domain);
euler = parcluster('local');
pool = parpool(euler,18);
T = 3*2*pi;
deltaT = T/15;
timeInterval = [0, T];
stepSize = deltaT;
range= 0:stepSize:T;
sens = zeros(length(initialPosition), length(range));

parfor i = 1:length(initialPosition)
   x0 = initialPosition(i,:);
   lDerivative = @(t,y) d_phi(t,y,0, false);
   [t, uncertAtTime] = computeAlongTrajectory(lDerivative, x0, timeInterval, stepSize );
   sens(i,:) = uncertAtTime;
   %lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
   %criteps2(i) = getCriticalEpsilon(lDerivative,x0,uncertAtTime, t, epsilons);
end

pool.delete();
save('SensitivityAlongTrajs_.mat', 'sens');
