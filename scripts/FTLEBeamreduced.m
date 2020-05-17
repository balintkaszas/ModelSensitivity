addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);

odegrad = generateFirstOrderODEGrad(beamStruct);
dygrad = @(t,x) odegrad(t,x,0,false);

initialPosition = zeros(32,1);
timeSpan = [0,0.001];
r = 5;
dT = 1e-8;
[cgmotd, ~] = computeOTD_separategradSerial(dy, dygrad, initialPosition.', timeSpan,r, dT);
[cgmfull, ~] = computeCGInvariants(dy, initialPosition.', timeSpan,'finiteDifference', false);
