addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x,e,eov) odefunction(t,x,e,eov);
timoshenkoBeam = DynSystem(dy, ndof_spv, [1,1]);
%%first run:
pool = parpool('local', 16);
mc = monteCarloUQ(timoshenkoBeam, 1e-5, 1024, true, 0);
save('mc1.mat', 'mc');
pool.delete();