addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');
odefunction = generateFirstOrderODE(beamStruct);


dy = @(t,x) TIMOSHENKODE(t,x.');

x0 = zeros(32,1);
endtime = 6.28;
dT = 0.001;
time = 0:dT:endtime;
%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;



nDeltas = 14;
nPoints = 5;
FTLES = zeros(nDeltas,nPoints);
Times = zeros(nDeltas,nPoints);

delt = linspace(-3,-16,nDeltas);
deltas = 10.^delt;


tol = 1e-13;
for i=1:nDeltas
    disp(i);
  
    for j = 1:nPoints
        delta = deltas(i);
        maxPlace = rand(32, 1)*1e-6;
        tic;
        [eigmax, ~] = computeCGInvariantsode15s_variabledelta(dy, maxPlace.', [0, 1.5], 'finiteDifference', false, tol, delta);
        Times(i,j) = toc;
        FTLES(i,j) = log(eigmax)/(3);
    end
end

save('ftle_timoshenko_sensitivity_to_delta_withode15s_file.mat', 'Times', 'FTLES','tol', 'deltas');

