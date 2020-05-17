
addPath;
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


rhs = A*X.' + Fnl;
%dy = @(t,y) odfunc(t, y, A, Fnn);


% compute the gradient symbolically 
syms t;
Fnn = matlabFunction(rhs, 'vars', {t,X}, 'File', 'TIMOSHENKODE');
%matlabFunction(dy, 'vars', {X}, );


X = [x;xd].';
MassMtx = [eye(ndof), zeros(ndof); zeros(ndof), M];
Alin = [zeros(ndof),eye(ndof);-K,-C];
Fnl = [zeros(ndof,1);-f];
syms t;
rhs = Alin*X.' + Fnl;
Fnn = matlabFunction(rhs, 'vars', {t,X}, 'File', 'TIMOSHENKODEWithMass');
MassMtx = double(MassMtx);
save('MMtx.mat', 'MassMtx');