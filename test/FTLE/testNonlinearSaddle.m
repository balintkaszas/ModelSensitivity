%%Testing the Implementation on a transformed nonlinear saddle problem. 
%%Setting up the parameters and initial conditions
domain = [-2, 2; 0, 4];
resolution = [100,100];
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
%[eigmaxFiniteDiff, traceFiniteDiff] = computeCGInvariants(derivative, initialPosition,  timespan, 'finiteDifference');
%ftl = reshape(log(eigmaxFiniteDiff)/(4), resolution);

%figure(3);
%surf(A1, A2,ftl);shading interp; axis equal;axis tight;colorbar; title('FTLE'); xlabel('u'); ylabel('v')
%view([0 0 1]); axis equal; axis tight; shading interp;camlight

% 
% 
% %%Calculating the model sensitivity
 %[MS, MSfull, hist] = modelSensitivityGlobal(derivative, initialPosition, resolution, timespan, t/20, 'finiteDifference');
 %dlmwrite('FTLE_eovSine_alt.txt', ftl, 'delimiter','\t','precision',16) %% Use every significant digit
 %dlmwrite('MS_FDSine_alt.txt', MS, 'delimiter','\t','precision',16) 
 %dlmwrite('MSFull_FDSine_alt.txt', MSfull, 'delimiter','\t','precision',16) 
 maxPlace = [1, pi/2];
 [t, uncertestimate] = computeAlongTrajectoryVariationGeneralFullmodel(derivative, maxPlace, [0,2], 2/20. );
 dlmwrite('time.txt', t.', 'delimiter','\t','precision',16) %% Use every significant digit

 dlmwrite('uncert_at_ridge.txt', uncertestimate.', 'delimiter','\t','precision',16) %% Use every significant digit
% 
% 
% figure(2);
% surf(A1, A2,MS);shading interp; axis equal;axis tight;colorbar; title('MS'); xlabel('u'); ylabel('v')
% view([0 0 1]); axis equal; axis tight; shading interp;camlight
% 
% figure(3);
% surf(A1, A2,MSfull);shading interp; axis equal;axis tight;colorbar; title('MS_\sigma'); xlabel('u'); ylabel('v')
% view([0 0 1]); axis equal; axis tight; shading interp;camlight
