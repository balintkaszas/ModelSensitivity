%%% Generate Duffing Model Sens. Fields
epsilon = 0.01;
omegapert = 10;
%%CDV params:
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;


dyref = @(t,x) d_charneyDeVore(t,x,p, false);
dygrad = @(t,x) d_charneyDeVore_grad(t,x,p);
%dypert = @(t,x) d_phi(t,x,0, false) + [0;epsilon*cos(omegapert*t)];
cdv = DynSystem(dyref, 6, [1,1], dygrad);

timeSpan1 = [0, 15];
k = rand(6,1)+1;
b0 = rand(6,1)+1;
b0 = b0./norm(b0);
maxPlace = [0.7864,0.,0.8848, 0, 0., 0.];

endtime = 15;
dT = 1e-4;
time = 0:dT:endtime;
[~,yfref] = ode45(dyref, time, maxPlace, odeset('RelTol', 1e-8));  %reference trajectory
epsilon = 1e-2;
[t,ms1] = modelSensitivityTrajectory(cdv, maxPlace, timeSpan1, 0.01, false, 0);
disp('MS1 Done');

sigma = eye(6)./6;
F = @(t,X) dyref(t,X) + epsilon*b0*sin(dot(k,X))*cos(omegapert*t);
G = @(t,X) epsilon*sigma;
obj = sde(F, G, 'StartState', maxPlace);
%rng('default');

[Paths,Times,Z] = simulate(obj,length(time)-1, 'DeltaTime' , dT, 'nTrials',2000);
diffs1 = Paths-yfref;

TimestoWrite = Times(1:100:end);
diffsToWrite = diffs1(1:100:end,:,:);
MSE1 = mean(sum(diffs1.^2,2),3);



save('CDVMC.mat', 't', 'TimestoWrite', 'MSE1','epsilon', 'ms1','diffsToWrite', 'maxPlace', 'b0', 'k', 'epsilon');