addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);
timoshenkoBeam = DynSystem(dy, ndof_spv, [1,1]);
%%first run:
pool = parpool('local', 16);
mc = monteCarloUQ(timoshenkoBeam, 1e-5, 2048, true, 3);
save('mcTRUE3.mat', 'mc');
pool.delete();