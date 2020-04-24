function [cgEigenvalueMax, cgTrace] = computeOTD_separategradSerial(derivative, derivativeEov, initialPosition, timeSpan,r, dT)
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

Tp = OTDEquation(derivative, derivativeEov, initialPosition, timeSpan, r, dT);

%loop through the array of nSystem by nSystem matrices and compute
%eigenvalues for each entry:
reducedEig = arrayfun(@(idx)svd(Tp(:,:,idx)),1:nRows,'UniformOutput',false); 
reducedEig = cell2mat(reducedEig);
cgEigenvalueMax = (reducedEig(1,:).^2); %%svd already returns singular values sorted
cgTraceReduced = sum(reducedEig.^2, 1); %sum over columns, gives the trace for the reduced operator
cgEigenvalueMin = reducedEig(end,:).^2; %to 'fill' in the remaining n-r dimensions, use the smallest eigenvalue
cgTrace = cgTraceReduced + (nSystem - r)*cgEigenvalueMin;
end


function Tp = OTDEquation(derivative, derivativeEov, initialPosition, timeSpan, r, dT)
    %% Calculates the Cauchy Green strain tensor from the equation of variations
    % append the initial nSystem x nRows to a nSystem + nSystem^2 x nRows
    % matrix

    nSystem =  size(initialPosition,2); % number of variables
    nRows = size(initialPosition,1); %number of gridpoints
    t0 = timeSpan(1);

    %% Setting up initial conditions for all fields
    U0 = zeros(nSystem, r, nRows); 
    tp0 = ones(1,r);
    kmax = (timeSpan(2) - timeSpan(1))/dT;
    Tp0 = eye(r);
    Tp0 = repmat(Tp0, 1,1, nRows);
    t = t0:dT:(kmax-1)*dT+t0;

    %dFlowMap0 = eye(nSystem);
    %dFlowMap0 = reshape(dFlowMap0,nSystem*nSystem,1); %% IC for EOV is unity
    %initGrid = [initialPosition,repmat(transpose(dFlowMap0),size(initialPosition,1),1)]; % prepare the ensemble ICs
    L = gradientL(derivativeEov, t0, initialPosition, nSystem, nRows); %% derivative at the initial time instant

    for i=1:nRows %% Setting up initial condition for the Reduced order dynamics
        L_sym = ( L(:,:, i) + L(:, :, i)' ) /2; %% symmetric part of the jacobian
        [V,Diag] = eig(L_sym);
        [~,IX] = sort(diag(Diag),'descend');
        U0(:,:, i) = V(:,IX(1:r));
    end
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
    for k=1:kmax %Runge-Kutta steps
        %% ------------- Direct Calculation of FTLE ------------------
        K1 = Velocity_Field(derivative, t(k), Xp, nSystem, nRows);
        K2 = Velocity_Field(derivative, t(k)+dT/2, Xp+dT*K1/2,  nSystem, nRows);
        K3 = Velocity_Field(derivative, t(k)+dT/2, Xp+dT*K2/2,  nSystem, nRows);
        K4 = Velocity_Field(derivative, t(k)+dT, Xp+dT*K3, nSystem, nRows);

        [K1_OTD, k1_otd] = Velocity_Field_OTD(derivativeEov, t(k), Xp,            U,               Tp, nSystem, nRows, r);
        [K2_OTD, k2_otd] = Velocity_Field_OTD(derivativeEov, t(k)+dT/2, Xp+dT*K1/2,   U+dT*K1_OTD/2,   Tp+dT*k1_otd/2, nSystem, nRows, r);
        [K3_OTD, k3_otd] = Velocity_Field_OTD(derivativeEov,   t(k)+dT/2, Xp+dT*K2/2,   U+dT*K2_OTD/2,   Tp+dT*k2_otd/2, nSystem, nRows, r);
        [K4_OTD, k4_otd] = Velocity_Field_OTD(derivativeEov, t(k)+dT, Xp+dT*K3,     U+dT*K3_OTD,     Tp+dT*k3_otd ,nSystem, nRows, r);

        Xp = Xp + dT*(K1+2*K2+2*K3+K4)/6;

        U  = U  + dT*(K1_OTD+2*K2_OTD+2*K3_OTD+K4_OTD)/6;

        Tp = Tp + dT*(k1_otd+2*k2_otd+2*k3_otd+k4_otd)/6;
    %% ----------------- Performing Gram-Schmidt -------------
        for p=1:nRows
         for i=1:r
            for j=1:i-1
                U(:,i,p) = U(:,i,p) - (U(:,i,p)'*U(:,j,p))*U(:,j,p);
            end
            U(:,i,p) = U(:,i,p)/sqrt(U(:,i,p)'*U(:,i,p));
         end
        end

    end
end


function grad = gradientL(f, t, initPositions, nSystem, nRows)
    %% f is expected to return nSystem by nSystem mtx
    % initPositions is nRows by nSystem 
    grad = zeros(nSystem, nSystem, nRows);
    for i = 1:nRows
        rowvalue = f(t, initPositions(i,:).');
        grad(:, :, i) = rowvalue;
    end
end


function dy = Velocity_Field(f, t, X, nSystem, nRows)
% function to compute the velocity vector for all points of the grid
    %% expects nSystem by nRows vector

dy = zeros(nRows, nSystem);
    for i=1:nRows
        fullState =  f(t, X(i,:).');
        dy(i, :) = fullState;
    end
end


function [R, R_otd]= Velocity_Field_OTD(f, t, X, U,Tp, nSystem, nRows, r)
    %% expects nSystem by nRows vector

L = gradientL(f, t, X, nSystem, nRows);
R = zeros(nSystem, r, nRows);
R_otd = zeros( r,r, nRows);
    for i=1:nRows
        R(:, :, i) = L(:, :,  i)*U( :, :, i);
        Q = eye(nSystem) - U(:, :, i)*U(:, :, i)';
        R(:, :, i) = Q*R(:, :, i);

        Lr = U(:, :, i)'*L(:, :, i)*U(:, :, i);
        temp = Lr*Tp(:, :, i);
        R_otd(:, :, i) = temp;
    end
end
