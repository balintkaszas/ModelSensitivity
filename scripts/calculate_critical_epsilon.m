lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
resolution = [100,100];
initialPosition = initialize_ic_grid(resolution,domain);

if isempty(gcp('nocreate'))
    parpool('local', 36);
end
dq = parallel.pool.DataQueue;
wb = waitbar(0, 'Pleasex wait...');
N = prod(resolution);
wb.UserData = [0 N];
afterEach(dq, @(varargin) iIncrementWaitbar(wb));
afterEach(dq, @(idx) fprintf('Completed iteration: %d\n', idx));

 
epsilons = logspace(-4, 2, 100);
criteps = zeros(length(initialPosition), 1);
T = 2*pi;
deltaT = T/100;
timeInterval = [0, T];
stepSize = deltaT;
length(initialPosition)

parfor i = 1:length(initialPosition)
   x0 = initialPosition(i,:);
   lDerivative = @(t,y) d_phi(t,y,0, false);
   [t, uncertAtTime] = computeAlongTrajectory(lDerivative, x0, timeInterval, stepSize );
   lDerivative = @(t,y,e, eov) d_phi(t,y,e, eov);
   criteps(i) = getCriticalEpsilon(lDerivative,x0,uncertAtTime, t, epsilons);
      send(dq, i);

end

save("criticalEps__g1.mat", "criteps");
close(wb);
pool.delete();

function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end

