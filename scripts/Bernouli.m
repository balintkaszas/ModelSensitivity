addPath;

ndof = 16;  % Even number for the beam
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

epsilon = 1e-6;
f = @(t,x) bernoullibeam(t, x, 0, false, A, Fnn, Fphi, ne);
g = @(t,x) epsilon;%*Fphi./Fphi;
disp(epsilon.*Fphi./norm(Fphi));
maxPlace = zeros(ndof_spv, 1);
Nsamples = 5000;
dt = 1e-6*pi;
time = 0:dt:pi;
[~,yfref] = ode45(@(t,y) bernoullibeam(t, y, 0, false, A, Fnn, Fphi, ne), [0,pi], maxPlace); 
plot(yfref(:, ndof-1), yfref(:, 2*ndof-1), '.')
xlabel('Position of the last node');
ylabel('Velocity of the last node');
BB = load('Bernoulli_MC_3_bigeps.mat');
Montecarlos = BB.Montecarlos;
ndof = BB.ndof;
hold on;
rng(0, 'twister'); %%for random samples
for i = 1:2
    r = randi([1 5000],1);
    yf = Montecarlos(:, :, r);
    plot(yf(:, ndof-1), yf(:, 2*ndof-1), '.-');
end