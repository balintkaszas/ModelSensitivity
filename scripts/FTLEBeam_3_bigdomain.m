addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

beamStruct = load('beam.mat');

odefunction = generateFirstOrderODEPert(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);

timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0]);
domain = [-10,10; -10, 10];
resolution = [20,20];
timeSpan = [0,3];
vars = [16, 32];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);
tolerance = 1e-12;

[eigmax, ~] = computeCGInvariants(dy, init.points, [0, 3],  'finiteDifference', true, 1e-10);
eigmax = reshape(eigmax, resolution);
save('mn___.mat', 'eigmax')

save('FTLETimoshenko_GOOOD_big_3__________.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
