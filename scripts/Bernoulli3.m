cla;clf;
addPath;
ndof = 3;  % Even number for the beam
ndof_spv = 2*ndof;

ne = ndof/2;
[M,C,K,fphi,fnl]=L_Bernoulli_Beam_Model(ne);

% ----- First order form -----
% \dot{x} = Ax + Fnl(x) + \epsilon*Fphi
A = [zeros(ndof),eye(ndof);-M\K,-M\C];
Fnl = [zeros(ndof,1);-M\fnl];
Fphi = [zeros(ndof,1);M\fphi];
%[t, sol] = ode45(@(t,x) bernoullibeam(t, x, 0, false, A, Fnl, Fphi), [0,10], zeros(20,1));
resolution = [251,251];

domain = [-1.5, 1.5;-1.5, 1.5]*1e-3;

initialPosition = initialize_ic_grid(resolution, domain, 2);
initgrid = zeros(length(initialPosition), ndof_spv);
initgrid(:, ndof) = initialPosition(:, 1);
initgrid(:, 2*ndof) = initialPosition(:, 2);
pool = parpool('local', 32);
lDerivative = @(t,x,e, eov) bernoullibeam(t, x, e, eov, A, Fnl, Fphi);
%[eigmax,~] = computeCGInvariants(lDerivative, initialPosition,[0,2*pi], 'eoV');
%ftle = log(eigmax)/(2*(2*pi));
MSFull = modelSensitivityGlobal(lDerivative, initgrid, resolution, [0,4*pi],0.2*pi, 'finiteDifference');
save('Bernoullibeam_3.mat', 'domain', 'ndof', 'MSFull', 'resolution');




