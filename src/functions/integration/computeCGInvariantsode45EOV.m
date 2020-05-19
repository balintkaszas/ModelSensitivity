function [cgEigenvalueMax, cgTrace] = computeCGInvariantsode45EOV(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance)
%% Computes the largest eigenvalue of the Cauchy Green tensor, along with the trace
% Arguments:
% derivative: function handle that returns RHS of the dynamical system 
% xdot = f(x, t, epsilon)
% Must be of the form: 
% derivative = @(t,x,e, useEoV) deriv(t,x,epsilon,useEoV, params);
%
% initialPosition : Grid of initial positions. 
%   number of rows: number of gridpoints in the computational domain
%   number of columns: number of equations
%
% timeSpan: the time interval of integration timeSpan = [t0, t1]
% 
% method: 'finiteDifference' - Evaluate the flow map's derivative using an
% auxiliary grid of relative spacing 1e-8.
% method: 'eoV' - equation of variations. Must use useEoV flag in
% derivative, which then has to return a n+n^2 dimensional vector of
% derivatives.
% 
% isParallel: if true, will compute the invariants over the grid in
% parallel.
% 
% Usage:
% [lambdaMax,Traces] =
% computeCGInvariants(derivative,initialPosition,timeSpan,'eoV', isParallel);
%[lambdaMax,Traces] =
% computeCGInvariants(derivative,initialPosition,timeSpan,'finiteDifference', isParallel);

    nSystem =  size(initialPosition,2); % number of variables
    nRows = size(initialPosition,1); %number of gridpoints

    deriv = @(t,x) derivative(t,x);
    dFlowmap = CGfromeov(deriv, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance);
    %loop through the array of nSystem by nSystem matrices and compute
    %eigenvalues for each entry:
    cgStrainD = arrayfun(@(idx)svd(dFlowmap(:,:,idx)),1:nRows,'UniformOutput',false); 
    cgStrainD = cell2mat(cgStrainD);
    cgEigenvalueMax = (cgStrainD(1,:).^2); %%svd already returns singular values sorted
    cgTrace = sum(cgStrainD.^2, 1); %sum over columns, gives the trace
end




function dflowmap = CGfromeov(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance)
       %% Eigenvectors from main grid
        nSystem = size(initialPosition,2);
        nRows = size(initialPosition,1);
        dflowmap = zeros(nSystem, nSystem, nRows);
        deov = @(t,x) derivativeTogether(t,x, derivative, derivativeEOV, nSystem);
        if isParallel == true
            dq = parallel.pool.DataQueue;

            wb = waitbar(0, 'Please wait...');
            % Use the waitbar's UserData to track progress
            wb.UserData = [0 nRows];
            afterEach(dq, @(varargin) iIncrementWaitbar(wb));
            afterEach(dq, @(idx) fprintf('Completed iteration: %d\n', idx));
            parfor i= 1:nRows %regular loop through ICS

                ic = initialPosition(i,:);
                initEov = eye(nSystem);
                x0 = [ic.'; initEov(:)];
                [~, sol]  = ode45(deov, [timeSpan(1), timeSpan(1) + 0.5*(timeSpan(2)-timeSpan(1)), timeSpan(2)], x0, odeset('relTol', tolerance));
                sol = sol(end,:);
                dflowmap(:,:,i) = reshape(sol(nSystem+1:end), [nSystem, nSystem]);
                send(dq, i);
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
    xvariable = x(1:nSystem);
    dyx = deriv(t, xvariable);
    Grad = derivgrad(t, xvariable);
    xmtx = reshape(x(nSystem+1:end), [nSystem, nSystem]);
    dymtx = Grad*xmtx;
    dy = [dyx; dymtx(:)];
end

function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end