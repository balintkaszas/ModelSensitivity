function [t, MSsqtraj] = modelSensitivityTrajectory(system, point, timeSpan, timeStep)
%Wrapper to compute the square of the Model Sensitivity
%Integrates the maximal eigenvalue and trace of the CG strain tensor along
%a trajectory, and returns the time dependence.
%   system is a DynSystem object. Has methods .rhs() and possiply
%   .gradrhs()
%   point is the Initial point for the trajectory, along which we calculate
%   the model sensitivity
%   timespan is [t0,t] : the time interval of the computation

derivative = @(t,x,e,eov) system.rhs(t,x, e, eov);  % have to keep this form for compute
[t, MSsqtraj] = computeMSAlongTrajectoryGeneral(derivative, point, timeSpan, timeStep, 'finiteDifference', system.deltas);
end

