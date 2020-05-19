addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');

timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0], dygrad);
domain = [-50, 50; -100, 100];
resolution = [100,100];
timeSpan = [0,3];
vars = [28, 32];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);
tolerance = 1e-12;
[ftle, ~] = FTLE(timoshenkoBeam, init, timeSpan, true, 1e-10, false, 0);

save('ModelSensTimoshenko_2_DOF_28_32.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
