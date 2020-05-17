
addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);

odegrad = generateFirstOrderODEGrad(beamStruct);
dygrad = @(t,x) odegrad(t,x,0,false);

initialPosition = zeros(32,1);
r = 5;
dT = 1e-8;

times = linspace(0.5*pi, 100*pi, 200);
times = times(1:10);
ftle = zeros(size(times));
ftleanal = zeros(size(times));

for i=1:length(times)
    disp(i);
    timeSpan = [0,times(i)];

    [cgmotd, ~] = computeOTD_separategradSerial(dy, dygrad, initialPosition.', timeSpan,r, dT);
    ftle(i) = log(cgmotd)/(2*times(i));
end
%hold on;
%plot(times, ftleanal, 'o-');

%plot(times, ftle, 'o-');
save('ftle_zero_timoshenko_OTD.mat', 'times', 'ftle');




