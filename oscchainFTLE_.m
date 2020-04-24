%%TEsting oscillator chain, against duffing oscillator.
%Duffing's parameter setup, considerable driving
addPath;
n = 3;
m = 1;
k = -1;
c = 0.15;
g = 0.5;
T = 1;
alpha = 0.3;
pool = parpool('local', 32);


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,1) = initialPosition(:,1);
initGrid(:,1+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof1_kneg_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('1-1');



resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,2) = initialPosition(:,1);
initGrid(:,2+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof2_kneg_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('1-2');


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,3) = initialPosition(:,1);
initGrid(:,3+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof2_kneg_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('1-3');


n = 3;
m = 1;
k = 1;
c = 0.15;
g = 0.5;
T = 1;
alpha = 0.3;


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,1) = initialPosition(:,1);
initGrid(:,1+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_1dof_kpos_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('2-1');


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,2) = initialPosition(:,1);
initGrid(:,2+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_2dof_kpos_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('2-2');



resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,3) = initialPosition(:,1);
initGrid(:,3+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_3dof_kpos_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('2-3');



n = 3;
m = 1;
k = 1;
c = 1;
g = 0.5;
T = 1;
alpha = 1e-5;


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,1) = initialPosition(:,1);
initGrid(:,1+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof1_SHparams_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('3-1');


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,2) = initialPosition(:,1);
initGrid(:,2+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof2_SHparams_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('3-2');


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,3) = initialPosition(:,1);
initGrid(:,3+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof3_SHparams_3.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('3-3');


n = 5;
m = 1;
k = 1;
c = 1;
g = 0.5;
T = 1;
alpha = 0.3;


resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,1) = initialPosition(:,1);
initGrid(:,1+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof1_kpos_c1_bigdriving_5.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('4-1');



resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,2) = initialPosition(:,1);
initGrid(:,2+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof2_kpos_c1_bigdriving_5.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('4-2');



resolution = [250, 250];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 2*n);

initGrid(:,3) = initialPosition(:,1);
initGrid(:,3+n) = initialPosition(:,2);


deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
[eigmax, ~] = computeCGInvariantsPar(deriv, initGrid, [0, 2*pi], 'finiteDifference');
save('FTLEfield_oscchain_dof3_kpos_c1_bigdriving_5.mat', 'eigmax', 'n', 'm', 'k', 'c', 'g', 'T', 'alpha');
disp('4-3');

pool.delete();
