addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');

timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0]);
domain = [-10,10; -10, 10];
resolution = [250,250];
timeSpan = [0,8];
vars = [ndof, 2*ndof];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);
tolerance = 1e-12;
ftle = FTLE(timoshenkoBeam, init, timeSpan, true, tolerance);

save('FTLETimoshenko_GOOOD_big_8.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
