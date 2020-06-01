function MSsq = modelSensitivity(system, grid, timeSpan, isParallel, useReduced, r)
%Wrapper to compute the square of the Model Sensitivity
%Integrates the maximal eigenvalue and trace of the CG strain tensor along
%trajectories. 
%   system is a DynSystem object. Has methods .rhs() and possiply
%   .gradrhs()
%   grid is a Grid object, grid.points is the coordinates
%   timespan is [t0,t] : the time interval of the computation
%   isParallel is a flag. if true, will use parallelized trajectory
%   advection
derivative = @(t,x) system.rhs(t,x);  % have to keep this form for compute
derivativeEov = @(t,x) system.gradRhs(t,x);
MSsq = computeModelSensitivity(derivative, derivativeEov, grid.points, timeSpan, 1e-6, 1, system.deltas, isParallel, useReduced, r);
MSsq = reshape(MSsq, grid.resolution);
end

