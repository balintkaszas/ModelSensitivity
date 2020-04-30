abc = DynSystem(@(t,x)d_abc(t,x,0,false), 3, [1,1]);
init = Grid(3, [1,2], [100,100], [0,2*pi; 0,2*pi], 1e-5);

[t, ms] = modelSensitivityTrajectory(abc, [0,0,0], [0,8], 0.1);
