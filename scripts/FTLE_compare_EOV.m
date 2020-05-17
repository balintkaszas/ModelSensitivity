addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');


dyMASS = @(t,x) TIMOSHENKODEWithMass(t,x.');
l = load('MMtx.mat');
MassMtx = l.MassMtx;

domain = [-10,10; -10, 10];
resolution = [10,10];
timeSpan = [0,3];
vars = [ndof, 2*ndof];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
icGrid = init.points;
pool = parpool('local', 32);
tolerance = 1e-12;
[ft1, ~] = computeCGInvariantsode15s(dy, icGrid, timeSpan, 'finiteDifference', false, 1e-12);
ft1 = reshape(ft1, resolution);

save('FTLETimoshenko_ode15_only.mat',  'ft1', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
