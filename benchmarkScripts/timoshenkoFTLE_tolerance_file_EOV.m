addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');



nTolerances = 1;
nPoints = 5;
FTLES = zeros(nTolerances,nPoints);
Times = zeros(nTolerances,nPoints);

maxPlace = zeros(32, 1);
maxPlace(16) = 0;
maxPlace(32) = 2;

for j = 1:nPoints
    disp(j);
    x0 = maxPlace + rand(32,1)*1e-5;
    tic;
    [eigmax, ~] = computeCGInvariantsode45EOV(dy, dygrad, x0.', [0, 3], false, 1e-12);
    Times(1,j) = toc;
    FTLES(1,j) = log(eigmax);
end
save('ftle_timoshenko_sensitivity_to_tolerance_withode45s_file_EOVV2_reltol_mplace_8.mat', 'Times', 'FTLES', 'maxPlace');



dy = @(t,x) TIMOSHENKODEWithMass(t,x.');
l = load('MMtx.mat');
MassMtx = l.MassMtx;

nTolerances = 14;
nPoints = 5;
FTLES = zeros(nTolerances,nPoints);
Times = zeros(nTolerances,nPoints);

for j = 1:nPoints
    disp(j);

    x0 = maxPlace + rand(32,1)*1e-5;
    tic;
    [eigmax, ~] = computeCGInvariantsode15s_parallel_massmtx(dy, MassMtx, x0.', [0, 3],  'finiteDifference', false, 1e-13);
    Times(1,j) = toc;
    FTLES(1,j) = log(eigmax);
end



save('ftle_timoshenko_sensitivity_to_tolerance_withode15s_par_file_FD_mplace_8.mat', 'Times', 'FTLES', 'maxPlace');

