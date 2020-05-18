addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE_Pert(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');

timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0], dygrad);
domain = [-3,3; -3, 3];
resolution = [100,100];
timeSpan = [0,4];
vars = [ndof, 2*ndof];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 36);
tolerance = 1e-8;
[ftle, ~] = FTLE(timoshenkoBeam, init, timeSpan, true, tolerance, false, 0);

save('FTLETimoshenko_1_DOF_pert.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
