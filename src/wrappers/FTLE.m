function ftle = FTLE(system, grid, timeSpan, isParallel, tolerance, method)
%FTLE wrapper to compute the finite time lyapunov exponent using a
%DynSystem
%object.
%   system is a DynSystem object. Has methods .rhs() and possiply
%   .gradrhs()
% grid is a Grid object, grid.points is the coordinates, grid.difference
%   is the distance between gridpoints.
% timespan is [t0,t] : the time interval of the computation
% isParallel is a flag. if true, will use parallelized trajectory
%   advection
% tolerance: absolute error tolerance for the integrator

    switch method
        case 'eov'

            derivative = @(t,x) system.rhs(t,x);  % have to keep this form for computeCGInvariants
            derivativeEOV = @(t,x) system.gradRhs(t, x);
            [cgEigmax, ~] = computeCGInvariantsEOV(derivative, derivativeEOV, grid.points, timeSpan, isParallel, tolerance); %keep only maximal eigenvalue from the outputs
        case 'finitedifference'
            derivative = @(t,x) system.rhs(t,x);  
            [cgEigmax, ~] = computeCGInvariantsFiniteDiff(derivative, grid.points, timeSpan, isParallel, grid.difference, tolerance); 
        otherwise
            disp('Unknown FTLE calculation method');
    end
    ftle = reshape(log(cgEigmax), grid.resolution)/(timeSpan(1)-timeSpan(0));
end

