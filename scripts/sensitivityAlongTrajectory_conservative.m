addpath(genpath('integration'));
addpath('models');
addpath('test');
lDerivative = @(t,y,e) d_phi_conservative(t,y,e);
resolution = [200,200];
domain = [-1.5, 1.5; -1.5, 1.5];

initialPosition = initialize_ic_grid(resolution,domain);

if isempty(gcp('nocreate'))
    parpool('local', 36);
end
dq = parallel.pool.DataQueue;
wb = waitbar(0, 'Please wait...');
N = prod(resolution);
wb.UserData = [0 N];
afterEach(dq, @(varargin) iIncrementWaitbar(wb));
afterEach(dq, @(idx) fprintf('Completed iteration: %d\n', idx));


T = 3*2*pi;
deltaT = T/100;
timeInterval = [0, T];
stepSize = deltaT;
range= 0:stepSize:T;
sens = zeros(length(initialPosition), length(range));

parfor i = 1:length(initialPosition)
   x0 = initialPosition(i,:);
   lDerivative = @(t,y) d_phi(t,y,0, false);
   [t, uncertAtTime] = computeAlongTrajectory(lDerivative, x0, timeInterval, stepSize );
   sens(i,:) = uncertAtTime;
    send(dq, i);
   %lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
   %criteps2(i) = getCriticalEpsilon(lDerivative,x0,uncertAtTime, t, epsilons);
end
close(wb);

save('SensitivityAlongTrajs____conservative__6pi_300.mat', 'sens');
function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end