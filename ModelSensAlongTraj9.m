addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE_Pert(t,x);
dygrad = @(t,x) TIMOSHENKODE_Grad_NL(t,x);

timoshenkoBeam = DynSystem(dy, ndof_spv, [1,1], dygrad);
pt = zeros(32,1);
pt(16) = 0.01;
pt(32) = -0.03;
timeSpan = [0,1];
timeStep = 0.05;

[t, ms2] = modelSensitivityTrajectory(timoshenkoBeam, pt.', timeSpan, timeStep, false, 0);


%plot(t, ms2, '--')
save('MSfineFirstHalf_16_32.mat','t', 'ms2', 'timeSpan');