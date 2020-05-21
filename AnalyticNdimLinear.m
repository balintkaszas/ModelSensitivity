%%% General case
addPath;
Amtx = [-2.3876,0.7226 ,2.1248; 0.2926 , -1.1591  , -0.9315; 2.4248  ,-0.9315 , -4.9887 ];
b = ones(3,1);
Ainvb = Amtx\b;
sigma = [0, 0, 0.8660; 0, 0.8660, 0; 0.8660, 0, 0.8660];




f0 = @(x) funzero(x, Amtx);
x0=eye(3);

endtime =20;
dT = 0.01;
time = 0:dT:endtime;
%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-6;

lm = eig(Amtx);
lmax = real(lm(3));

sigmanorm = trace(sigma.'*sigma);
Symma = (Amtx.' + Amtx);
traceSymmaInv = trace(inv(Symma));
MSasympt = epsilon^2*norm(b).^2./norm(lmax)^2 - epsilon^2*sigmanorm*traceSymmaInv;
perturbvector = zeros(32,1);
perturbvector(16:end) = 1;
F = @(t,X) Amtx*X + epsilon*b;
G = @(t,X) epsilon*sigma;
x0 = zeros(3,1);
obj = sde(F, G, 'StartState', x0);
rng('default');

[Paths,Times,Z] = simulate(obj,2000, 'DeltaTime' , dT, 'nTrials', 1000);
MSE = squeeze(mean(sum(Paths.^2,2),3));
save('PathsNEW.mat', 'Paths', 'Times', 'MSE', 'MSasympt');



b2 = norm(b)*real(V(:,3))/norm(real(V(:,3)));


x0 = zeros(3,1);
F = @(t,X) Amtx*X + epsilon*b2;
G = @(t,X) epsilon*sigma;
obj = sde(F, G, 'StartState', x0);
rng('default');

[Paths,Times,Z] = simulate(obj,2000, 'DeltaTime' , dT, 'nTrials', 1000);
MSE = squeeze(mean(sum(Paths.^2,2),3));
save('PathsNEWOpt.mat', 'Paths', 'Times', 'MSE', 'MSasympt');



b2 = norm(b)*real(V(:,3))/norm(real(V(:,3)));


x0 = zeros(3,1);
F = @(t,X) Amtx*X+epsilon*b2;
Times = 0:dT:20;
[~, sol] = ode45(F, Times, x0, odeset('absTol', 1e-12));
MSE = squeeze(sum(sol.^2,2));
MSasympt = epsilon^2*norm(b).^2./norm(lmax)^2;
save('PathsNEWOptDet.mat', 'Times', 'MSE', 'MSasympt');



