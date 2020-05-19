addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;

dy = @(t,x) TIMOSHENKODE(t,x.');
dygrad = @(t,x) TIMOSHENKODE_GRAD(t,x.');


dyMASS = @(t,x) TIMOSHENKODEWithMass(t,x.');
l = load('MMtx.mat');
MassMtx = l.MassMtx;


oneIC = rand(32,1);
twoIC = rand(32,1);

t = linspace(0, 10, 1000);
[~, sol45] = ode45(dy, t, oneIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, sol15] = ode15s(dy, t, oneIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, solmass] = ode15s(dyMASS, t, oneIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Mass', MassMtx));
disp('oneRef')

oneICPert = oneIC + 1e-5*rand(32,1);
[~, sol45PERT] = ode45(dy, t, oneICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, sol15PERT] = ode15s(dy, t, oneICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, solmassPERT] = ode15s(dyMASS, t, oneICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Mass', MassMtx));
disp('onePert')

t = linspace(0, 10, 1000);
[~, sol452] = ode45(dy, t, twoIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, sol152] = ode15s(dy, t, twoIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, solmass2] = ode15s(dyMASS, t, twoIC, odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Mass', MassMtx));
disp('twoRef')

twoICPert = twoIC + 1e-5*rand(32,1);
[~, sol452PERT] = ode45(dy, t, twoICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, sol152PERT] = ode15s(dy, t, twoICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12));
[~, solmass2PERT] = ode15s(dyMASS, t, twoICPert, odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'Mass', MassMtx));
disp('twoPert')


save('SOLUTIONS_3_solvers.mat', 't', 'oneIC', 'twoIC', 'sol45', 'sol15', 'solmass', 'sol45PERT', 'sol15PERT', 'solmassPERT', 'sol452', 'sol152', 'solmass2', 'sol452PERT', 'sol152PERT', 'solmass2PERT');