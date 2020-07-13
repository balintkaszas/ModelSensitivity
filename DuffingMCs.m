%%% Generate Duffing Model Sens. Fields

omegapert = 10;
dyref = @(t,x) d_phi(t,x,0, false);
dygrad = @(t,x) d_phi_grad(t,x);
%dypert = @(t,x) d_phi(t,x,0, false) + [0;epsilon*cos(omegapert*t)];
duff = DynSystem(dyref, 2, [1,1], dygrad);

timeSpan1 = [0, 6*pi];


p1= [-0.8253; -0.48795];
p2 = [-1.0904;-0.87349];
p3 = [0.70482; -0.24699];

endtime = 6*pi;
dT = 1e-4*pi;
time = 0:dT:endtime;
[~,yfref] = ode45(dyref, time, p1, odeset('RelTol', 1e-8));  %reference trajectory
epsilon = 1e-2;
[t,ms1] = modelSensitivityTrajectory(duff, p1, timeSpan1, 0.01*pi, false, 0);
disp('MS1 Done');


F = @(t,X) dyref(t,X) + [0;epsilon*cos(omegapert*t)];
G = @(t,X) [0;epsilon];
obj = sde(F, G, 'StartState', p1);
%rng('default');

[Paths,Times,Z] = simulate(obj,length(time)-1, 'DeltaTime' , dT, 'nTrials',2000);
diffs1 = Paths-yfref;
MSE1 = mean(sum(diffs1.^2,2),3);





[~,yfref] = ode45(dyref, time, p2, odeset('RelTol', 1e-8));  %reference trajectory
epsilon = 1e-2;
[t,ms2] = modelSensitivityTrajectory(duff, p2, timeSpan1, 0.01*pi, false, 0);
disp('MS2 Done');


F = @(t,X) dyref(t,X) + [0;epsilon*cos(omegapert*t)];
G = @(t,X) [0;epsilon];
obj = sde(F, G, 'StartState', p2);
%rng('default');

[Paths,Times,Z] = simulate(obj,length(time)-1, 'DeltaTime' , dT, 'nTrials',2000);
diffs2 = Paths-yfref;
MSE2 = mean(sum(diffs2.^2,2),3);


[~,yfref] = ode45(dyref, time, p3, odeset('RelTol', 1e-8));  %reference trajectory
epsilon = 1e-2;
[t,ms3] = modelSensitivityTrajectory(duff, p3, timeSpan1, 0.01*pi, false, 0);
disp('MS3 Done');


F = @(t,X) dyref(t,X) + [0;epsilon*cos(omegapert*t)];
G = @(t,X) [0;epsilon];
obj = sde(F, G, 'StartState', p3);
%rng('default');

[Paths,Times,Z] = simulate(obj,length(time)-1, 'DeltaTime' , dT, 'nTrials',2000);
diffs3 = Paths-yfref;
MSE3 = mean(sum(diffs3.^2,2),3);

save('DuffingMC.mat', 't', 'time', 'MSE1', 'MSE2', 'MSE3', 'epsilon', 'ms1', 'ms2', 'ms3', 'diffs1', 'diffs2', 'diffs3', 'p1', 'p2', 'p3');