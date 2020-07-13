function MS = computeModelSensitivity(derivative, derivgrad, initialPosition, timeSpan, toleranceFD, toleranceInt, Deltas, isParallel, method)
%% Computes the model sensitivity, by integrating the invatiants of the CG tensor, _along_ trajectories, started from an initial grid.

% toleranceFD: absolute error tolerance for the Runge-Kutta integrators
% toleranceInt: Error tolerance for the computation of the integrals in the
%   integral() method
% 
% Deltas: Deltas = [D1, D2]. Best guess on the magnitudes of the perturbations. Must be a 2 by
% 1 vector. D1: Magnitude of the deterministic perturbation
%   D2: !Square! of the magnitude of the stochastic perturbation.
% Either one can be zero.
%
% method: 'eov' or 'finitedifference'
% Usage:
% 
% returns an array that has the resolution provided as input, with each entry corresponding to the (local) model sensitivity

t0 = timeSpan(1);
t = timeSpan(2);
nRows = size(initialPosition, 1);
nSystem = size(initialPosition, 2);
Integral = zeros(nRows, 2);
if isParallel == true
    parfor i = 1:nRows
        ic = initialPosition(i,:);
        funsq = @(s) vecFtlePullback(derivative, derivgrad, ic, s, t, t0, toleranceFD, method);
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
MS = (D1^2)*Integral(:,1).^2 + D2*Integral(:,2);  %% this returns the model sensitivity!
end




function [sqrtEig,cgTrace] = ftlePullback(dy, derivgrad,  x0, s, t, t0, toleranceFD, method)
    % Computes the invariants of the Cauchy-Green strain tensor from s to
    % t, at a point x(s), on a trajectory starting from x0 at t0.
    % Not parallel. 
    %
    % Advects the IC from t0 to s
    % calculates the CG tensor's invariants on the interval [s, t]
    % two special cases:    s == t0, no need to advect
    %                       s = t, the flowmap is the identity tensor
    %                       cgMax = 1, cgTrace = dim.
    x = x0;
    if(s ~= t0)
        x = ode45Ensemble(dy, [t0, s], x0, false); % Advect the gridpoints to time s
    end
    if(s ~= t)
        % here, we must set isPalallel to false.
        [cgMax, cgTrace] = computeCGInvariantsWrapper(dy, derivgrad, x, [s, t], false, toleranceFD, 1e-3, method);
    else
        cgMax = 1;
        cgTrace = numel(x0);  %% dimension of the system
    end
    sqrtEig = sqrt(cgMax);  % computeCGInvariants returns the eigenvalue, need to integrate the square root of it

end

function resultVec = vecFtlePullback(dy, derivgrad,  x0, s, t, t0, toleranceFD, method)
    % vectorize the function ftlePullback in the variable s.
    % has the same arguments, accepts an array as "s".
    % returns a 2 by 1 array
    nMax = length(s);
    [sqrteig, cgtrace] = arrayfun(@(idx)ftlePullback(dy, derivgrad, x0, s(idx), t, t0, toleranceFD, method),1:nMax); 
    resultVec(1) = sqrteig;
    resultVec(2) = cgtrace;
end
