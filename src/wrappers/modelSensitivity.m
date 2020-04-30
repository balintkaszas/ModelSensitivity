function MSsq = modelSensitivity(system, grid, timeSpan, timeStep, isParallel)
%Wrapper to compute the square of the Model Sensitivity
%Integrates the maximal eigenvalue and trace of the CG strain tensor along
%trajectories. 
%   system is a DynSystem object. Has methods .rhs() and possiply
%   .gradrhs()
%   grid is a Grid object, grid.points is the coordinates
%   timespan is [t0,t] : the time interval of the computation
%   isParallel is a flag. if true, will use parallelized trajectory
%   advection
derivative = @(t,x,e,eov) system.rhs(t,x, e, eov);  % have to keep this form for compute
MSsq = modelSensitivityGlobal(derivative, grid.points, grid.resolution,  timeSpan, timeStep, 'finiteDifference', isParallel, system.deltas); %keep only maximal eigenvalue from the outputs
end

