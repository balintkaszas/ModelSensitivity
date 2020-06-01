addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE_Pert(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');

timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0], dygrad);
domain = [-500, 500; -500, 500; -500, 500];
resolution = [50,50,50];
timeSpan = [0,12];
vars = [16,23, 32];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);
tolerance = 1e-12;
[ftle, ~] = FTLE(timoshenkoBeam, init, timeSpan, true, 1e-11, false, 0);

save('FTLETimoshenko_2_DOF_16_23_32_resonance.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
