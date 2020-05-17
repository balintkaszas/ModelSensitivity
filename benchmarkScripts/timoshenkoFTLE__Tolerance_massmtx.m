
addPath;
ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;


dy = @(t,x) TIMOSHENKODEWithMass(t,x.');
l = load('MMtx.mat');
MassMtx = l.MassMtx;
x0 = zeros(32,1);


nTolerances = 11;
nPoints = 5;
FTLES = zeros(nTolerances,nPoints);
Times = zeros(nTolerances,nPoints);

tolerances = linspace(-6,-16,11);
tolerances = 10.^tolerances;



for i=1:nTolerances
    disp(i);
    for j = 1:nPoints
        tol = tolerances(i);
        maxPlace = rand(32, 1)*1e-6;
        tic;
        [eigmax, ~] = computeCGInvariantsode15s_parallel_massmtx(dy, MassMtx, maxPlace.', [0, 1.5], 'finiteDifference', false, tol);
        Times(i,j) = toc;
        FTLES(i,j) = log(eigmax)/(3);
    end
end

save('ftle_timoshenko_sensitivity_to_tolerance_withMassMtx.mat', 'Times', 'FTLES');

