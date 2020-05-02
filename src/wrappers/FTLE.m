function ftle = FTLE(system, grid, timeSpan, isParallel)
%FTLE wrapper to compute the finite time lyapunov exponent using a
%DynSystem
%object.
%   system is a DynSystem object. Has methods .rhs() and possiply
%   .gradrhs()
%   grid is a Grid object, grid.points is the coordinates
%   timespan is [t0,t] : the time interval of the computation
%   isParallel is a flag. if true, will use parallelized trajectory
%   advection
derivative = @(t,x) system.rhs(t,x);  % have to keep this form for computeCGInvariants
[Eigmax, ~] = computeCGInvariants(derivative, grid.points, timeSpan, 'finiteDifference', isParallel); %keep only maximal eigenvalue from the outputs
ftle = reshape(log(Eigmax), grid.resolution)/(2*(timeSpan(2) - timeSpan(1)));
end

