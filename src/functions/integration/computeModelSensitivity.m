function MS = computeModelSensitivity(derivative, derivgrad, initialPosition, timeSpan, toleranceFD, toleranceInt, Deltas, isParallel, useReduced, r)
%% Computes the model sensitivity, by integrating the invatiants of the CG tensor, _along_ trajectories, started from an initial grid.
% creates a subdivision of the time interval, and evaluates the invariants
% over time intervals of increasing length [t_0, tMax], [t_1, tMax], ...

% Arguments:
% derivative: function handle that returns RHS of the dynamical system 
% xdot = f(x, t, epsilon)
% Must be of the form: 
% derivative = @(t,x,e, useEoV) deriv(t,x,epsilon,useEoV, params);
%
% initialPosition : Grid of initial positions. 
%   number of rows: number of gridpoints in the computational domain
%   number of columns: number of equations
% resolution: resolution of the initial grid. Needed for the final
% reshaping
%
% timeSpan: the time interval of integration timeSpan = [t0, t1]
% stepSize: intermediate steps to take
% method: 'finiteDifference' - Evaluate the flow map's derivative using an
% auxiliary grid of relative spacing 1e-8.
% method: 'eoV' - equation of variations. must use useEoV flag in
% derivative, which then has to return a n+n^2 dimensional vector of
% derivatives.
%
% isParallel: if true, will compute the invariants over the grid in
% parallel.
%
% Deltas: Deltas = [D1, D2]. Best guess on the magnitudes of the perturbations. Must be a 2 by
% 1 vector. D1: Magnitude of the deterministic perturbation
% D2: !Square! of the magnitude of the stochastic perturbation.
% Either one can be zero.
%
%
% Usage:
% MS =
% modelSensitivityGlobal(derivative, initialPosition, resolution, timeSpan, stepSize, 'eoV', isParallel, [D1, D2])
% returns an array that has the resolution provided as input, with each entry corresponding to the (local) model sensitivity
% MS is the model sensitivity in case of purely deterministic error
t0 = timeSpan(1);
t = timeSpan(2);
nRows = size(initialPosition, 1);
nSystem = size(initialPosition, 2);
Integral = zeros(nRows, 2);
if isParallel == true
    parfor i = 1:nRows
        ic = initialPosition(i,:);
        funsq = @(s) vecFtlePullback(derivative, derivgrad, ic, s, t, t0, toleranceFD, useReduced, r);
        Integral(i, :) =  integral(funsq, t0, t, 'ArrayValued', true, 'AbsTol', toleranceInt);
    end
else
    for i = 1:nRows
        ic = initialPosition(i,:);
        funsq = @(s) vecFtlePullback(derivative, derivgrad, ic, s, t, t0, toleranceFD, useReduced, r);
        Integral(i, :) =  integral(funsq, t0, t, 'ArrayValued', true, 'AbsTol', toleranceInt);
    end
end


D1 = Deltas(1); %get weights from the argument
D2 = Deltas(2);
% return the value in the correct shape
MS = (D1^2)*Integral(:,1).^2 + D2*Integral(:,2);  %% this returns the square of the model sensitivity!
end




function [sqrteig,cgtrace] = ftlePullback(dy, derivgrad,  x0, s, t, t0, toleranceFD, useReduced, r)
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
        x = ode45_vector_nonvectorized(dy, [t0, s], x0, false); % Advect the gridpoints to time s
    end
    if(s ~= t)
        if useReduced == true
            dT = (t-s)/10;
            [cgmax, cgtrace] = computeOTD_ODE45(dy, derivgrad, x, [s, t], r, dT, true, false);
        else
            [cgmax, cgtrace] = computeCGInvariantsode45EOV(dy, derivgrad, x, [s, t], false, toleranceFD);
        end
    else
        cgmax = 1;
        cgtrace = numel(x0); 
    end
    sqrteig = sqrt(cgmax);  % computeCGInvariants returns the eigenvalue, need to integrate the square root of it

end

function resultVec = vecFtlePullback(dy, derivgrad,  x0, s, t, t0, toleranceFD, useReduced, r)
    % vectorize the function ftlePullback in the variable s.
    % has the same arguments, accepts an array as "s".
    % returns a 2 by 1 array
    nMax = length(s);
    [sqrteig, cgtrace] = arrayfun(@(idx)ftlePullback(dy, derivgrad, x0, s(idx), t, t0, toleranceFD, useReduced, r),1:nMax); 
    resultVec(1) = sqrteig;
    resultVec(2) = cgtrace;
end




