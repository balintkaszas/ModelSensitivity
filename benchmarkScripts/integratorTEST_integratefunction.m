addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;
knull = zeros(6,6);

dy = @(t,x) d_charneyDeVore(t,x,p, false, 0,knull,0);

%f = @(t,x) d_charneyDeVore(t,x,p, false, epsilon,k1, .1);

nTimesteps = 10;
timesteps = linspace(10, 1000, nTimesteps);
nPoints = 5;
MS = zeros(nTimesteps,nPoints);
Times = zeros(nTimesteps,nPoints);
timesteps = 1./timesteps;
t0 = 0;
t = 1.5;
rng(0);
x0 = rand(6,1)*1e-1;
%asd2
funsq = @(s) vecFtlePullbackFOR(dy, x0.', s, t, t0);
funsqAppl = @(s) vecFtlePullback(dy, x0.', s, t, t0);
disp('\n First')
tic;
%s = 0:0.1:1;
%funsq(s)
Integral =  integral(funsq, t0, t, 'ArrayValued', true, 'AbsTol', 1e-5);
toc

disp('\n First')
tic;
%s = 0:0.1:1;
%funsq(s)
Integral2 =  integral(funsqAppl, t0, t, 'ArrayValued', true, 'AbsTol', 1e-5);
toc

tic;
dt = 0.01;
old = modelSensitivityGlobal(dy, x0.', [1,1], [t0, t], dt, 'finiteDifference', false, [1,1]);
toc
%Integral
disp(Integral(1).^2 + Integral(2))
disp(Integral2(1).^2 + Integral2(2))

disp(old)




%save('MS_duffing_sensitivity_to_timestep.mat', 'Times', 'MS', 'timesteps');

function [sqrteig,cgtrace] = ftlePullback(dy, x0, s, t, t0)
    % Computes the invariants of the Cauchy-Green strain tensor from s to
    % t, at a point x(s), on a trajectory starting from x0 at t0.
    % Advects the IC from t0 to s
    % calculates the CG tensor's invariants on the interval [s, t]
    % two special cases:s == t0, no need to advect
    % s = t, the flowmap is the identity tensor
    % cgmax = 1, cgtrace = dim.
    % for now, it only works for one point. 
    x = x0;

    if(s ~= t0)
        x = ode45_vector_nonvectorized(dy, [t0, s], x0, 'false'); % Advect the gridpoints to time s
    end
    if(s ~= t)
        [cgmax, cgtrace] = computeCGInvariantsode15s(dy, x, [s, t], 'finiteDifference', 'false', 1e-12);
    else
        cgmax = 1;
        cgtrace = numel(x0); 
    end
    sqrteig = sqrt(cgmax);  % computeCGInvariants returns the eigenvalue, need to integrate the square root of it

end

function resultVec = vecFtlePullback(dy, x0, s, t, t0)
    % vectorize the function ftlePullback in the variable s.
    % has the same arguments, accepts an array as "s".
    % returns a 2 by 1 array
    nMax = length(s);
    [sqrteig,cgtrace]  = arrayfun(@(idx)ftlePullback(dy, x0, s(idx), t, t0),1:nMax, 'UniformOutput', false); 
    sqrteig = cell2mat(sqrteig);
    cgtrace = cell2mat(cgtrace);
    resultVec(1) = sqrteig;
    resultVec(2) = cgtrace;
end


function resultVec = vecFtlePullbackFOR(dy, x0, s, t, t0)
    % vectorize the function ftlePullback in the variable s.
    % has the same arguments, accepts an array as "s".
    % returns a 2 by 1 array
    nMax = length(s);
    resultVec = zeros(2,nMax);
    for i = 1:nMax
        [resultVec(1,i), resultVec(2,i)] = ftlePullback(dy, x0, s(i), t, t0);
    end

end

