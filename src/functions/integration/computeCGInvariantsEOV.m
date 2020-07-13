function [cgEigenvalueMax, cgTrace] = computeCGInvariantsEOV(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance)
%% Computes the largest eigenvalue of the Cauchy Green tensor, along with the trace from the equation of variations
% Arguments:
% derivative: function handle that returns RHS of the dynamical system 
%   xdot = f(x, t)
%   Must be of the form: 
%   derivative = @(t,x) deriv(t,x, params);
% derivativeEOV: function handle to compute Jacobian of RHS of the
%   dynamical system. Must return a matrix 
% 
% initialPosition : Grid of initial positions. 
%   number of rows: number of gridpoints in the computational domain
%   number of columns: number of equations
%
% timeSpan: the time interval of integration timeSpan = [t0, t1]
% 
% isParallel: if true, will compute the invariants over the grid in
% parallel.
% 
% tolerance: absolute error tolerance for the integrator
% 
% Usage:
% [lambdaMax,Traces] =
% computeCGInvariantsEOV(derivative,derivativeEOV, initialPosition, timeSpan,isParallel);

    nSystem =  size(initialPosition,2); % number of variables
    nRows = size(initialPosition,1); %number of gridpoints

    deriv = @(t,x) derivative(t,x);
    dFlowmap = CGfromeov(deriv, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance);
    %loop through the array of nSystem by nSystem matrices and compute
    %eigenvalues for each entry:
    cgSvd = arrayfun(@(idx)svd(dFlowmap(:,:,idx)),1:nRows,'UniformOutput',false); 
    cgSvd = cell2mat(cgSvd);
    cgEigenvalueMax = (cgSvd(1,:).^2); %%svd already returns singular values sorted
    cgTrace = sum(cgSvd.^2, 1); %sum over columns, gives the trace
end




function dflowmap = CGfromeov(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance)
       %% Eigenvectors from main grid
        nSystem = size(initialPosition,2);
        nRows = size(initialPosition,1);
        dflowmap = zeros(nSystem, nSystem, nRows);
        deov = @(t,x) derivativeTogether(t,x, derivative, derivativeEOV, nSystem);
        if isParallel == true

            parfor i= 1:nRows %parallel loop through ICS
                ic = initialPosition(i,:);
                initEov = eye(nSystem);
                x0 = [ic.'; initEov(:)];
                %use intermediate timestep to save memory.
                [~, sol]  = ode45(deov, [timeSpan(1), timeSpan(1) + 0.5*(timeSpan(2)-timeSpan(1)), timeSpan(2)], x0, odeset('relTol', tolerance));
                sol = sol(end,:);
                dflowmap(:,:,i) = reshape(sol(nSystem+1:end), [nSystem, nSystem]);
            end   
            
        else
            for i= 1:nRows %regular loop through ICS
                ic = initialPosition(i,:);
                initEov = eye(nSystem);
                x0 = [ic.'; initEov(:)];
                [~, sol]  = ode45(deov, [timeSpan(1), timeSpan(1) + 0.5*(timeSpan(2)-timeSpan(1)), timeSpan(2)], x0, odeset('relTol', tolerance));
                sol = sol(end,:);
                dflowmap(:,:,i) = reshape(sol(nSystem+1:end), [nSystem, nSystem]);
            end
        end
end

function dy = derivativeTogether(t,x, deriv, derivgrad, nSystem)
%Computes the derivatives of the full system and the linearized system
%simultaneously. 
%x: nSystem+nSystem*nSystem length vector
%deriv: rhs of the original system
%derivgrad: jacobian
    xvariable = x(1:nSystem);
    dyx = deriv(t, xvariable);
    Grad = derivgrad(t, xvariable);
    xmtx = reshape(x(nSystem+1:end), [nSystem, nSystem]);
    dymtx = Grad*xmtx;
    dy = [dyx; dymtx(:)]; %returns a vector of the same shape as x.
end

