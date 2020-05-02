function ftleR = FTLEReduced(system, grid, timeSpan, r, dT, isParallel)
%Wrapper to use the OTD (optimally time dependent) modes to compute the finite time lyapunov exponent using a
%DynSystem 
%object. It approximates the CG tensor's eigenvectors with a set of r, optimally
%selected vectors.
%   system is a DynSystem object. Has methods .rhs() and must have the
%   method .gradrhs() too.
%   grid is a Grid object, grid.points is the coordinates
%   timespan is [t0,t] : the time interval of the computation
%   r is the order of the OTD subspace reduction
%   dT is the timestep for the OTD computation 
%   isParallel is a flag. if true, will use parallelized trajectory
%   advection
derivative = @(t,x) system.rhs(t, x);  % have to keep this form for computeCGInvariants
derivativeGrad = @(t, x) system.gradrhs(t, x);

if isParallel == true
    [Eigmax, ~] = computeOTD_separategrad(derivative, derivativeGrad, grid.points, timeSpan, r, dT); %keep only maximal eigenvalue from the outputs
else
    [Eigmax, ~] = computeOTD_separategradSerial(derivative, derivativeGrad, grid.points, timeSpan, r, dT); %keep only maximal eigenvalue from the outputs
end

ftleR = reshape(log(Eigmax), grid.resolution)/(2*(timeSpan(2) - timeSpan(1)));
end

