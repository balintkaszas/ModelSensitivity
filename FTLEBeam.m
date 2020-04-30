addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x,e,eov) odefunction(t,x,e,eov);
timoshenkoBeam = DynSystem(dy, ndof_spv, [0,0]);
domain = [-4*1e-4, 4*1e-4; -4*1e-4, 4*1e4];
resolution = [250,250];
timeSpan = [0,10];
vars = [ndof, 2*ndof];
init = Grid(ndof_spv, vars, resolution, domain, 1e-8);
pool = parpool('local', 32);

ftle = FTLE(timoshenkoBeam, init, timeSpan, true);

save('FTLETimoshenko.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
pool.delete();