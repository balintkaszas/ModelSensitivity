addpath(genpath('integration'));
addpath('models');
addpath('test');
domain = [-1.5, 1.5; -1.5, 1.5];
resolution = [200,200];

if isempty(gcp('nocreate'))
    parpool('local', 36);
end
dq = parallel.pool.DataQueue;
wb = waitbar(0, 'Please wait...');
N = prod(resolution);
wb.UserData = [0 N];
afterEach(dq, @(varargin) iIncrementWaitbar(wb));
afterEach(dq, @(idx) fprintf('Completed iteration: %d\n', idx));


initialPosition = initialize_ic_grid(resolution,domain);
n = length(initialPosition);
criticalEpsilon = zeros(n,1);
devs = zeros(n, length(0:0.025:3));
addpath('models');
addpath('integration');
addpath('integration/grid');
parfor i = 1:n
    x0 = initialPosition(i,:);
    [devs(i,:), criticalEpsilon(i)] = criticalEps(@(t,x,e,eov) d_phi(t,x,e,eov),[0, 2*pi],x0);
    send(dq, i);
end
save('criticalEpsLinearity3.mat','criticalEpsilon');
save('EpsLinearitydependenceLong.mat','devs');



close(wb);
pool.delete();

function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end

