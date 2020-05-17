
addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');
odefunction = generateFirstOrderODE(beamStruct);

A = generateFirstOrderODEMTX(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);
timoshenkoBeam = DynSystem(dy, ndof_spv, [1,1]);
%%first run:

deriv = @(t,x) timoshenkoBeam.rhs(t,x);    %% unperturbed system
%maxPlace = rand(system.dimension, 1)*1e-8;
maxPlace = zeros(timoshenkoBeam.dimension, 1);

endtime = 6.28;
dT = 0.001;
time = 0:dT:endtime;
%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;



maxplaces = rand(32,100)*1e-6;
mss = zeros(100,1);
times = zeros(100,1);

for i=1:length(maxplaces)
    disp(i);
    maxplace = maxplaces(:,i);
    tic;
    ms = modelSensitivityGlobal(deriv, maxPlace.',[1,1], [0, 3], 0.1, 'finiteDifference', false, [1,1]);
    times(i) = toc;
    mss(i) = ms;
end

save('mss_s_timoshenko.mat', 'times', 'mss');

