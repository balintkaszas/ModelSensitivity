%%Testing the Implementation on a transformed nonlinear saddle problem. 
%%Setting up the parameters and initial conditions
domain = [-2, 2; -4, 4];
resolution = [100,200];
initialPosition = initialize_ic_grid(resolution,domain,2);
A1 = initialPosition(:,1);
A2 = initialPosition(:,2);

A1 = reshape(A1, resolution);
A2 = reshape(A2, resolution);

t0 = 0;
t = 2;


timespan = [0, 2];
derivative = @(t,x,e,eov)testSystemSine(t,x,e,eov);

%% Numerical FTLE field and trace. 
[eigmaxFiniteDiff, traceFiniteDiff] = computeCGInvariants(derivative, initialPosition,  timespan, 'finiteDifference');
ftl = reshape(log(eigmaxFiniteDiff)/(4), resolution);


%%Calculating the model sensitivity, from equation of variations
[MSeov, MSfulleov] = modelSensitivityGlobal(derivative, initialPosition, resolution, timespan, t/100, 'eoV');
dlmwrite('FTLE_finiteDiff_sine.txt', ftl, 'delimiter','\t','precision',16) %% Use every significant digit
dlmwrite('MS_eov_sine.txt', MS, 'delimiter','\t','precision',16) 
dlmwrite('MSFull_eov_sine.txt', MSfull, 'delimiter','\t','precision',16) 


%%Calculating the model sensitivity, from equation of variations
[MSfd, MSfullfd] = modelSensitivityGlobal(derivative, initialPosition, resolution, timespan, t/100, 'finiteDifference');
dlmwrite('MS_finiteDiff_sine.txt', MS, 'delimiter','\t','precision',16) 
dlmwrite('MSFull_finiteDiff_sine.txt', MSfull, 'delimiter','\t','precision',16) 

