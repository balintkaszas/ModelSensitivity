function [cgEigenvalueMax,cgTrace] = computeOTD_ODE45(derivative, derivativeEov, initialPosition, timeSpan,r, dT, withTrace, isParallel)
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
% method: 'eoV' - equation of variations. must use useEoV flag in
% derivative, which then has to return a n+n^2 dimensional vector of
% derivatives.
% method: 'advect' - calculate largest LE by the growth of random
% perturbations. Following the algorithm of Benettin et al., only computes
% the largest LE!
%
% Usage:
% [lambdaMax,Traces] =
% computeCGInvariants(derivative,initialPosition,timeSpan,'eoV');
%[lambdaMax,Traces] =
% computeCGInvariants(derivative,initialPosition,timeSpan,'finiteDifference');
%[lambdaMax,Traces] =
% computeCGInvariants(derivative,initialPosition,timeSpan,'advect');

nSystem =  size(initialPosition,2); % number of variables
nRows = size(initialPosition,1); %number of gridpoints
Tp = zeros(r,r, nRows);
Determinant = zeros(1, nRows);

if isParallel == true
    parfor i = 1:nRows
        ic = initialPosition(i,:);
        [Tp(:,:,i), Determinant(i)] = OTDEquation(derivative, derivativeEov, ic, timeSpan, r, dT, withTrace);
    end
end
if isParallel == false
    for i = 1:nRows
        ic = initialPosition(i,:);
        [Tp(:,:,i), Determinant(i)] = OTDEquation(derivative, derivativeEov, ic, timeSpan, r, dT, withTrace);
    end
end
%loop through the array of nSystem by nSystem matrices and compute
%eigenvalues for each entry:
reducedEig = arrayfun(@(idx)svd(Tp(:,:,idx)),1:nRows,'UniformOutput',false); 
reducedEig = cell2mat(reducedEig);
cgEigenvalueMax = (reducedEig(1,:).^2); %%svd already returns singular values sorted

cgProduct = prod(reducedEig.^2, 1); %Product of the first r eigenvalues
cgTraceReduced = sum(reducedEig.^2,1);
lambdaRPlus1 = (Determinant.^2./cgProduct).^(nSystem-r);%%Estimate for the first eigenvalue that was NOT computed
%to 'fill' in the remaining n-r dimensions, use the smallest eigenvalue
cgTrace = cgTraceReduced + (nSystem - r)*lambdaRPlus1;
end


function [Tp,Determinant] = OTDEquation(derivative, derivativeEov, initialPosition, timeSpan, r, dT, withTrace)
    %% Calculates the Cauchy Green strain tensor from the equation of variations
    % append the initial nSystem x nRows to a nSystem + nSystem^2 x nRows
    % matrix

    nSystem =  size(initialPosition,2); % number of variables
    t0 = timeSpan(1);

    %% Setting up initial conditions for all fields
    U0 = zeros(nSystem, r); 
    tp0 = ones(1,r);
    kmax = (timeSpan(2) - timeSpan(1))/dT;
    Tp0 = eye(r);
    
    %dFlowMap0 = eye(nSystem);
    %dFlowMap0 = reshape(dFlowMap0,nSystem*nSystem,1); %% IC for EOV is unity
    %initGrid = [initialPosition,repmat(transpose(dFlowMap0),size(initialPosition,1),1)]; % prepare the ensemble ICs
    L = derivativeEov(t0, initialPosition); %% derivative at the initial time instant

    L_sym =( L+ L.' ) /2; %% symmetric part of the jacobian
    [V,Diag] = eig(L_sym);
    [~,IX] = sort(diag(Diag),'descend');
    U0 = V(:, IX(1:r));
    %Setting up output fields:
    Xp = initialPosition;
    %disp(size(initGrid));
    %XpwEov = initGrid;
    U  = U0;
    Tp = Tp0;
    % Because of the added condition of orthogonality at each step, we
    % perform the RK iteration explicitly.
    %Evolve the trajectories, the basis and the reduced basis
    %simultaneously
    Determinant = 1;
    for k=1:kmax %Runge-Kutta steps
        %% ------------- Direct Calculation of FTLE ------------------
        currtime = (k-1)*dT;
        [~, solXp] = ode45(derivative, [currtime, currtime+dT/2, currtime+dT], Xp, odeset('absTol', 1e-8));
        Xp = solXp(end, :);
        icU = [U(:); Tp(:)];
        if withTrace == true
            icU = [U(:); Tp(:); Determinant];
        end
        derivv = @(t,x) Velocity_Field_OTD(derivativeEov, t, Xp, x, nSystem, r, withTrace);
        [~, solUT] = ode45(derivv, [currtime, currtime+dT/2, currtime+dT], icU, odeset('absTol', 1e-8));
        
        UTnew = solUT(end, :);
        Unew = UTnew(1:numel(U));
        Tpnew = UTnew(numel(U)+1:end);
        U = reshape(Unew, size(U));
        Tp = reshape(Tpnew, size(Tp));
        Determinant = 0;
        if withTrace == true
            Determinant = solUT(end);
        end
    %% ----------------- Performing Gram-Schmidt -------------
     for i=1:r
        for j=1:i-1
            U(:,i) = U(:,i) - (U(:,i)'*U(:,j))*U(:,j);
        end
        U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
     end

    end
end




function output = Velocity_Field_OTD(derivEov, t, X, xin, nSystem, r, withTrace)
    %% expects nSystem by nRows vector
Uflat = xin(1:nSystem*r);
U = reshape(Uflat, [nSystem, r]);

Tflat = xin(nSystem*r+1:end);
if withTrace == true
    Tflat = xin(nSystem*r+1:end-1);
    Determinant = xin(end);
end
Tp = reshape(Tflat, [r,r]);
L = derivEov(t, X);
R = zeros(nSystem, r);
R_otd = zeros(r,r);
R = L*U;
Q = eye(nSystem) - U*U';
R = Q*R;

Lr = U'*L*U;
temp = Lr*Tp;
R_otd = temp;
output = [R(:);R_otd(:)];
if withTrace == true
    output = [R(:);R_otd(:);trace(L)*Determinant];
end
end