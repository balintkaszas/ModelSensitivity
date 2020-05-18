addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

[odefunction, odegrad, A] = generateFirstOrderODE(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);
dygrad = @(t,x) odegrad(t,x,0,false);

point = zeros(32,1);%*1e-8;
timeSpan = [0,6];
dT = 1e-7;
r = 5;
asdd = DynSystem(dy, 32, [1,1], dygrad);

%time = linspace(0, 20, 100);
%p2 = point;
%p2(16) = 1e-10;
%[~,sol1] = ode15s(dy2, time, point, odeset('relTol', 1e-13)); 
%[~,sol2] = ode15s(dy2, time, p2, odeset('relTol', 1e-13)); 
%diff = sol1 - sol2;

%plot(sol1 - sol(2))
%eig(A)
[t2, u2] = modelSensitivityTrajectory(asdd, point.', timeSpan, 3/10);
%[Eigmax, ~] = computeOTD_separategradSerial(dy, dygrad, point.', timeSpan, r, dT); %keep only maximal eigenvalue from the outputs
%ftleR = log(Eigmax)/(2*(timeSpan(2) - timeSpan(1)));

%[Eigmaxo, ~] = computeCGInvariants(dy, point.', timeSpan, 'finiteDifference', false); %keep only maximal eigenvalue from the outputs
%ftle = log(Eigmaxo)/(2*(timeSpan(2) - timeSpan(1)));

%disp('done');

