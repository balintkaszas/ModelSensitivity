function MS = modelSensitivityGlobal2(derivative, initialPosition, resolution, timeSpan, stepSize, method, isParallel, Deltas)
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



nSystem =  size(initialPosition,2); % number of variables
nRows = size(initialPosition,1); %number of gridpoints

grid = initialPosition;
range = (timeSpan(1):stepSize:timeSpan(2)-stepSize); %% specify the subdivision of the interval
n = length(range);
CGMaxToSum = zeros(nRows, size(range,1)); %% predefine the matrix containting CG maximal eigenvalues, for each subinterval and each gridpoint
CGTraceToSum = zeros(nRows, size(range,1)); %% same for the CG trace
for i = 1:n %%loop over the time interval

    [cgmax, cgtrace] = computeCGInvariants(derivative, grid, [range(i), timeSpan(2)], method, isParallel); %Just pass the method argument to computeCGInvariants
    grid = ode45_vector_nonvectorized(@(t,y)derivative(t,y,0,false), [range(i), range(i) + stepSize], grid, isParallel); % Advect the gridpoints to the next time instant 
    CGMaxToSum(:,i) = sqrt(cgmax); % computeCGInvariants returns the eigenvalue, need to integrate the square root of it
    CGTraceToSum(:,i) = cgtrace;  
end

summCGMax = trapz(range, CGMaxToSum, 2); %Trapezoid integration separately
summCGTrace = trapz(range, CGTraceToSum, 2); %then add them together with the weights

D1 = Deltas(1); %get weights from the argument
D2 = Deltas(2);
% return the value in the correct shape
MS = reshape((D1^2)*summCGMax.^2 + D2*summCGTrace, fliplr(resolution));  %% this returns the square of the model sensitivity!
end




