%% gradmtx

%addPath;

ndof = 16;  %16 degree of freedom Timoshenko beam
ndof_spv = 2*ndof;
system = load('beam.mat');

C = system.C;
K = system.K;
M = system.M;
f = system.f;
x = system.x;
xd = system.xd;
ndof = length(x);

% \dot{X} = AX + Fnl(X) 
% X = [x, xdot]
X = [x;xd];

A = [zeros(ndof),eye(ndof);-M\K,-M\C];
Fnl = [zeros(ndof,1);-M\f];
omega = 0.3;
Fphi = zeros(ndof_spv, 1);
Fphi(end) = 1;
syms t;

size(X)
rhs = A*X + Fnl + 1e-4*Fphi*sin(omega*t);
%dy = @(t,y) odfunc(t, y, A, Fnn);
