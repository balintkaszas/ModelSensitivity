function [cgEigenvalueMax, cgTrace] = computeCGInvariantsFiniteDiff(derivative, initialPosition, timeSpan, isParallel, difference, tolerance)
%% Computes the largest eigenvalue of the Cauchy Green tensor, along with the trace from a finite difference approximation
% Arguments:
% derivative: function handle that returns RHS of the dynamical system 
%   xdot = f(x, t)
%   Must be of the form: 
%   derivative = @(t,x) deriv(t,x, params);
%
% initialPosition : Grid of initial positions. 
%   number of rows: number of gridpoints in the computational domain
%   number of columns: number of equations
%
% timeSpan: the time interval of integration timeSpan = [t0, t]
% 
% isParallel: if true, will compute the invariants over the grid in
% parallel.
% 
% difference: distance between gridpoints
%
% tolerance: absolute error tolerance for the integrator
% 
% Usage:
% [lambdaMax,Traces] =
% computeCGInvariantsEOV(derivative,derivativeEOV, initialPosition, timeSpan, isParallel, tolerance);


    nSystem =  size(initialPosition,2); % number of variables
    nRows = size(initialPosition,1); %number of gridpoints
    relDelta = 1e-5; %% use 1e-5 times closer auxiliary gridpoints for finite differencing

    dFlowmap = CGfromFD(derivative, initialPosition, difference, relDelta, timeSpan, isParallel, tolerance);
    %loop through the array of nSystem by nSystem matrices and compute
    %eigenvalues for each entry:
    cgSvd = arrayfun(@(idx)svd(dFlowmap(:,:,idx)),1:nRows,'UniformOutput',false); 
    cgSvd = cell2mat(cgSvd);
    cgEigenvalueMax = (cgSvd(1,:).^2); %%svd already returns singular values sorted
    cgTrace = sum(cgSvd.^2, 1); %sum over columns, gives the trace
end




function dflowmap = CGfromFD(derivative, initialPosition, Diff, relDelta, timeSpan, isParallel, tolerance)
       %% Eigenvectors from main grid
        nSystem = size(initialPosition,2);
        nRows = size(initialPosition,1);
        Delta = Diff*relDelta; 
        dflowmap = zeros(nSystem, nSystem, nRows);
        if isParallel == true
            parfor i= 1:nRows %parallel loop through ICS
                ic = initialPosition(i,:);
                for j = 1:nSystem
                    dir = zeros(size(ic));
                    dir(j) = 1; % perturbation vector, with only the jth coordinate being nonzero
                    ic1 = ic + dir * Delta; % new auxiliary gridpoint in each direction
                    ic2 = ic - dir * Delta; % for central difference, each direction has 2 new points
                    [~,sol1] = ode45(derivative, timeSpan, ic1, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    [~,sol2] = ode45(derivative, timeSpan, ic2, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    sol1 = sol1(end, :);
                    sol2 = sol2(end, :);
                    for k = 1:nSystem
                        dflowmap(j,k,i) = (sol1(k) - sol2(k))./(2*Delta); %central difference
                    end
                end
            end 
            
        else
            for i= 1:nRows %regular loop through ICS
                ic = initialPosition(i,:);
                for j = 1:nSystem
                    dir = zeros(size(ic));
                    dir(j) = 1; % perturbation vector, with only the jth coordinate being nonzero
                    ic1 = ic + dir * Delta; % new auxiliary gridpoint in each direction
                     ic2 = ic - dir * Delta; 
                    [~,sol1] = ode45(derivative, timeSpan, ic1, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    [~,sol2] = ode45(derivative, timeSpan, ic2, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    sol1 = sol1(end, :);
                    sol2 = sol2(end, :);
                    for k = 1:nSystem
                        dflowmap(j,k,i) = (sol1(k) - sol2(k))./(2*Delta); %central difference
                    end
                end
            end
        end
end

