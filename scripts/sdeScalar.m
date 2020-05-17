%%SDE

addPath;
%[~,yfref] = ode45(deriv, time, maxPlace);  %reference trajectory
epsilon = 1e-3;


a = -2;
b = 1;
sigma = 1;
F = @(t,X) a*X + epsilon*b;
G = @(t,X) epsilon*sigma;
x0 = 0;
obj = sde(F, G, 'StartState', x0);
rng('default');
dy = @(t,x) linscalar(t,x,a);
scalarSystem = DynSystem(dy, 1, [1,1]);
[t, mssq] = modelSensitivityTrajectory(scalarSystem, [0], [0, 3], 1/100);

[Paths,Times,Z] = simulate(obj, 6*1e2, 'DeltaTime' , 1e-2, 'nTrials', 2000);
MSE = mean(Paths.^2,3);
analit = analytic(Times, epsilon,b, a, sigma);


save('linearscalar.mat', 'Paths', 'Times', 't', 'mssq', 'analit');
function anal = analytic(t, epsilon,b, a, sigma) %% analytic expression for ms error.
    det = (epsilon.^2*b.^2/a.^2)*(exp(a*t) - 1).^2;
    stoch = ((epsilon^2*sigma^2)/(2*a))*(exp(2*a*t)-1);
    anal = det + stoch;
end


function dy = linscalar(t,x, a)
    dy = a*x;
end