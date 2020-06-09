addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE_Pert(t,x);
dygrad = @(t,x) TIMOSHENKODE_Grad_NL(t,x);
omega0 = 11.125;
T0 = 2*pi/omega0;
timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0], dygrad);
domain = [-1000, 1000; -1000, 1000];
resolution = [50,50];
timeSpan = [0,4*T0];
vars = [16,32];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);
tolerance = 1e-12;
[ftle, ~] = FTLE(timoshenkoBeam, init, timeSpan, true, 1e-8, false, 0);

save('FTLETimoshenko_2dd_2_DOF_16_32__resonance__biig.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();
