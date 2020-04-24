addpath(genpath('integration'));
addpath('models');
addpath('test');
domain = [-1.5, 1.5; -1.5, 1.5];
resolution = [200,200];
initialPosition = initialize_ic_grid(resolution,domain);
n = length(initialPosition);


T = 2*pi;
deltaT = T/100;
timeInterval = [0, T];
stepSize = deltaT;
range= 0:stepSize:T;


maxdevs = zeros(n,length(range));
addpath('models');
addpath('integration');
addpath('integration/grid');
if isempty(gcp('nocreate'))
    parpool('local', 36);
end
dq = parallel.pool.DataQueue;
wb = waitbar(0, 'Pleasex wait...');
N = prod(resolution);
wb.UserData = [0 N];
afterEach(dq, @(varargin) iIncrementWaitbar(wb));
afterEach(dq, @(idx) fprintf('Completed iteration: %d\n', idx));


parfor i = 1:n
    x0 = initialPosition(i,:);
    maxdevs(i,:) = maxDev(@(t,x,e,false) d_phi(t,x,e,false),x0, range, 0.005);
    send(dq, i);
end

save('maxdevs______2pi_100_timedep_VERYSMALLepsilon.mat','maxdevs');
close(wb);
pool.delete();

function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end

