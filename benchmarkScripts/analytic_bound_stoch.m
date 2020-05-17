addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODEPert(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);
timoshenkoBeam = DynSystem(dy, ndof_spv, [1,1]);
%%first run:

deriv = @(t,x) timoshenkoBeam.rhs(t,x);    %% unperturbed system
%maxPlace = rand(system.dimension, 1)*1e-8;
maxPlace = zeros(timoshenkoBeam.dimension, 1);

endtime = 6.28;
dT = 0.001;
time = 0:dT:endtime;
[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;

[t, mssq] = modelSensitivityTrajectory(timoshenkoBeam, maxPlace.', [0,2*pi], pi/10.);
save('nonautonomous_zero.mat', 't', 'mssq');


perturbvector = zeros(32,1);
perturbvector(16:end) = 1;
F = @(t,X) deriv(t,X) + perturbvector*epsilon;
G = @(t,X) epsilon*perturbvector;
x0 = zeros(32,1);
obj = sde(F, G, 'StartState', x0);
rng('default');

[Paths,Times,Z] = simulate(obj, 6.28/0.001, 'DeltaTime' , 0.001, 'nTrials', 4000);
save('Montecarlos_autonomous.mat', 'Times', 'Paths', 'epsilon');