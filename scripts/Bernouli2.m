cla;clf;
addPath;
ndof = 4;  % Even number for the beam
ndof_spv = 2*ndof;

ne = ndof/2;
[M,C,K,fphi,fnl]=L_Bernoulli_Beam_Model(ne);

% ----- First order form -----
% \dot{x} = Ax + Fnl(x) + \epsilon*Fphi
A = [zeros(ndof),eye(ndof);-M\K,-M\C];
Fnl = [zeros(ndof,1);-M\fnl];
Fphi = [zeros(ndof,1);M\fphi];
var = symvar(Fnl);
Fnn = matlabFunction(Fnl, 'vars', var);
resolution = [100,100];
domain = [-1.5, 1.5;-1.5, 1.5]*1e-2;

initialPosition = initialize_ic_grid(resolution, domain, 2);
initgrid = zeros(length(initialPosition), ndof_spv);
initgrid(:, ndof) = initialPosition(:, 1);
initgrid(:, 2*ndof) = initialPosition(:, 2);
%parpool('local', 8);
lDerivative = @(t,x,e, eov) bernoullibeam(t, x, e, eov, A, Fnn, Fphi, ne);
[eigmax,~] = computeCGInvariants(lDerivative, initgrid,[0,2*pi], 'finiteDifference');
%ftle = log(eigmax)/(2*(2*pi));
%save('Bernoullibeam_FTLEGood.mat', 'domain', 'ndof', 'eigmax', 'resolution');

%MSFull = modelSensitivityGlobal(lDerivative, initgrid, resolution, [0,2*pi],0.2*pi, 'finiteDifference');
%save('Bernoullibeam_MS.mat', 'domain', 'ndof', 'MSFull',  'resolution');
