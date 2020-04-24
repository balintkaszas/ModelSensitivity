function MS = modelSensitivityGlobalOTD2(derivative, derivativeEov, initialPosition, resolution, timeSpan, stepSize, r, isParallel, Deltas)
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
% method: 'advect' - calculate largest LE by the growth of random
% perturbations. Following the algorithm of Benettin et al., only computes
% the largest LE!
%
% Usage:
% [MS, MSFull] =
% modelSensitivityGlobal(derivative, initialPosition, resolution, timeSpan, stepSize, 'eoV')
% returns an array that has the resolution provided as input, with each entry corresponding to the (local) model sensitivity
% MS is the model sensitivity in case of purely deterministic error
% MSFull is the model sensitivity in case of both deterministic and
% stochastic errors. 


nSystem =  size(initialPosition,2); % number of variables
nRows = size(initialPosition,1); %number of gridpoints

grid = initialPosition;
range = (timeSpan(1):stepSize:timeSpan(2)-stepSize); %% specify the subdivision of the interval
n = length(range);
CGMaxToSum = zeros(nRows, size(range,1)); %% predefine the matrix containting CG maximal eigenvalues, for each subinterval and each gridpoint
CGTraceToSum = zeros(nRows, size(range,1)); %% same for the CG trace
%[dq, wbk] = progressBar(n);
for i = 1:n %%loop over the time interval
    %send(dq, i);
    %disp(i)
    if isParallel == true
        [cgmax,cgtrace] = computeOTD_separategrad(derivative, derivativeEov, grid, [range(i), timeSpan(2)], r, stepSize/10.); %Just pass the method argument to computeCGInvariants
    else
        [cgmax,cgtrace] = computeOTD_separategradSerial(derivative, derivativeEov, grid, [range(i), timeSpan(2)], r, stepSize/10.); %Just pass the method argument to computeCGInvariants
    end
    grid = ode45_vector_nonvectorized(derivative, [range(i), range(i) + stepSize], grid, isParallel); % Advect the gridpoints to the next time instant 
    CGMaxToSum(:,i) = sqrt(cgmax); 
    CGTraceToSum(:,i) = cgtrace;  %invariants over a subinterval
end
%CGMaxToSum(:,n+1) = 1.; 
%CGTraceToSum(:,n+1) = nSystem;  %invariants over a subinterval

summCGMax = trapz(range, CGMaxToSum, 2);
summCGTrace = trapz(range, CGTraceToSum, 2); % integrate separately, then add them together 

D1 = Deltas(1); %get weights from the argument
D2 = Deltas(2);
% return values
%MS = reshape(summCGMax, fliplr(resolution));
MS = reshape((D1^2)*summCGMax.^2 + D2*summCGTrace,fliplr(resolution)); 
%history = reshape(CGMaxToSum, [resolution(1), resolution(2), n]); 
%close(wb);
end




