
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
X = [x;xd].';

A = [zeros(ndof),eye(ndof);-M\K,-M\C];
Fnl = [zeros(ndof,1);-M\f];
omega = 0.3;
Fphi = zeros(ndof_spv, 1);
Fphi(end) = 1;
syms t;
rhs = A*X.' + Fnl + 1e-4*Fphi*sin(omega*t);
%dy = @(t,y) odfunc(t, y, A, Fnn);


% compute the gradient symbolically 
Fnn = matlabFunction(rhs, 'vars', {t,X}, 'File', 'TIMOSHENKODE_Pert');
%matlabFunction(dy, 'vars', {t, X} );
A = double(A);
save('A.mat', 'A');

Fgrad = jacobian(Fnl, X);
Fnn = matlabFunction(Fgrad, 'vars', {X}, 'File', 'TIMOSHENKODE_Grad_NL');

% 
% X = [x;xd].';
% MassMtx = [eye(ndof), zeros(ndof); zeros(ndof), M];
% Alin = [zeros(ndof),eye(ndof);-K,-C];
% Fnl = [zeros(ndof,1);-f];
% syms t;
% rhs = Alin*X.' + Fnl;
% Fnn = matlabFunction(rhs, 'vars', {t,X}, 'File', 'TIMOSHENKODEWithMass');
% MassMtx = double(MassMtx);
% save('MMtx.mat', 'MassMtx');