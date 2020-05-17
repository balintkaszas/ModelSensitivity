%%SDE

addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');

maxPlace = rand(32, 1).'*1e-6;
opts = odeset('relTol', 1e-12, 'absTol', 1e-12, 'Stats','on');

%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;
disp('ODE 45');
tic;
[~,yfref] = ode45(dy, time, maxPlace, opts);  %reference trajectory
toc

disp('ODE 15s');
tic;
[~,yfref] = ode15s(dy, time, maxPlace, opts);  %reference trajectory
toc

disp('ODE 23s');
tic;
[~,yfref] = ode23s(dy, time, maxPlace, opts);  %reference trajectory
toc

disp('ODE 113');
tic;
[~,yfref] = ode113(dy, time, maxPlace, opts);  %reference trajectory
toc

