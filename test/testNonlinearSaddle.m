%%Testing the Implementation on a transformed nonlinear saddle problem, the Sine-Ridge 
%%Setting up the parameters and initial conditions
domain = [-2, 2; 0, 4];
resolution = [50,50];

%Specifying the derivative and the jacobian for the system
derivative = @(t,x)testSystemSine(t,x);
derivativegrad = @(t,x)testSystemSine_grad(t,x);

%set stochastic term to 0. 
TestSystem = DynSystem(derivative, 2, [1,1], derivativegrad);

%setting up the grid:
ics = Grid(2, [1,2], resolution, domain, 1e-3);

timeSpan = [0,2];
%% Numerical FTLE field, both finite diff and EOV calculations
ftle = FTLE(TestSystem, ics,  timeSpan, false, 1e-7, 'finitedifference');
dlmwrite('FTLE_finitediff.dat', ftle, 'delimiter','\t','precision',16) %% Use every significant digit
ftle = FTLE(TestSystem, ics,  timeSpan, false, 1e-7, 'eov');
dlmwrite('FTLE_eov.dat', ftle, 'delimiter','\t','precision',16) %% Use every significant digit
disp('FTLE done')
%% Numerical MS field
ms = modelSensitivity(TestSystem, ics,  timeSpan, false, 1e-7, 'finitedifference');
dlmwrite('MS_finitediff.dat', ms, 'delimiter','\t','precision',16) %% Use every significant digit
ms = modelSensitivity(TestSystem, ics,  timeSpan, false, 1e-7, 'eov');
dlmwrite('MS_eov.dat', ms, 'delimiter','\t','precision',16) %% Use every significant digit
