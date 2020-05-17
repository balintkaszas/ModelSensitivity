addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');

nTimesteps = 6;
timesteps = linspace(200, 300, nTimesteps);
nPoints = 3;
MS = zeros(nTimesteps,nPoints);
Times = zeros(nTimesteps,nPoints);
timesteps = 1./timesteps;


pool = parpool('local', 32);

for i=1:nTimesteps
    disp(i);
    for j = 1:nPoints
        dt = timesteps(i);
        maxPlace = rand(32, 1)*1e-5;
        tic;
        ms = modelSensitivityGlobal(dy, maxPlace.', [1,1], [0, 1], dt, 'finiteDifference', false, [1,1]);
        Times(i,j) = toc;
        MS(i,j) = ms;
    end
end

save('MS_timoshenko_sensitivity_to_timestepSMALL_parallel.mat', 'Times', 'MS', 'timesteps');
pool.delete();
