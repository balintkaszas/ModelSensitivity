function [cgEigenvalueMax, cgTrace] = computeCGInvariants(derivative, initialPosition, timeSpan, method, isParallel, tolerance)
%%We trust this function... Tested against analytic example.
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

switch method
    case 'finiteDifference'
        derivEov = @(t,x) derivative(t,x);
        relDelta = 1e-5; %% convenience
        differences = diff(initialPosition);
        maximalDifference = max(differences, [], 'all'); %% max. difference over the grid
        Diff = 1e-3; %% set a default value in this case

        dFlowmap = CGfromFD(derivEov, initialPosition, Diff, relDelta, timeSpan, isParallel, tolerance);
end
    %loop through the array of nSystem by nSystem matrices and compute
    %eigenvalues for each entry:
    cgStrainD = arrayfun(@(idx)svd(dFlowmap(:,:,idx)),1:nRows,'UniformOutput',false); 
    cgStrainD = cell2mat(cgStrainD);
    cgEigenvalueMax = (cgStrainD(1,:).^2); %%svd already returns singular values sorted
    cgTrace = sum(cgStrainD.^2, 1); %sum over columns, gives the trace
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
                    [~,sol1] = ode15s(derivative, timeSpan, ic1, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    [~,sol2] = ode15s(derivative, timeSpan, ic2, odeset('relTol', 1e-12, 'absTol', tolerance)); 
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
                    [~,sol1] = ode15s(derivative, timeSpan, ic1, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    [~,sol2] = ode15s(derivative, timeSpan, ic2, odeset('relTol', 1e-12, 'absTol', tolerance)); 
                    sol1 = sol1(end, :);
                    sol2 = sol2(end, :);
                    for k = 1:nSystem
                        dflowmap(j,k,i) = (sol1(k) - sol2(k))./(2*Delta); %central difference
                    end
                end
            end
        end
end

