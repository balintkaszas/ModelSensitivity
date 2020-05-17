%%SDE

addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
beamStruct = load('beam.mat');

odefunction = generateFirstOrderODE(beamStruct);
dy = @(t,x) odefunction(t,x,0,false);

maxPlace = rand(32, 1)*1e-6;
time = [0, 3];
opts = odeset('relTol', 1e-12, 'absTol', 1e-12, 'Stats','on');

%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;
disp('ODE 45');
tic;
[~,yfref] = ode45(dy, time, maxPlace, opts);  %reference trajectory
toc

disp('ODE 15s');
tic;
[~,yfref] = ode15s(dy, time, maxPlace.', opts);  %reference trajectory
toc






%%% With mass matrix

C = beamStruct.C;
K = beamStruct.K;
M = beamStruct.M;
f = beamStruct.f;
x = beamStruct.x;
xd = beamStruct.xd;
ndof = length(x);
X = [x;xd];
MassMtx = [eye(ndof), zeros(ndof); zeros(ndof), M];
Alin = [zeros(ndof),eye(ndof);-K,-C];
Fnl = [zeros(ndof,1);-f];
Fnn = matlabFunction(Fnl, 'vars', X);
Alin = double(Alin);
odf = @(t,y) odfunc(t,y, Alin, Fnn);
opts = odeset('relTol', 1e-12, 'absTol', 1e-12, 'Mass', double(MassMtx), 'Stats','on');

disp('ODE 15s');
tic;
[~,yfref] = ode15s(odf, time, maxPlace, opts);  %reference trajectory
toc
function odf = odfunc(t, xvar, A, Fnn)
    unpacked = num2cell(xvar);
    odf = A*xvar + Fnn(unpacked{:}); % Fnn expects separate variables as argument, instead of a vector..
end
