function odeFunction = generateFirstOrderODE(system)
%from the system class that contains symbolic functions describing the
%second order system, it returns the right hand side of the first order
%system

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
Fnn = matlabFunction(Fnl, 'vars', X);
A = double(A);
Asparse = sparse(A);
odeFunction = @(t, xvar, e, eov) odfunc(t, xvar, e, eov, Asparse, Fnn);
end

function odf = odfunc(t, xvar, e, eov, A, Fnn)
    unpacked = num2cell(xvar);
    odf = A*xvar + Fnn(unpacked{:}); % Fnn expects separate variables as argument, instead of a vector..
end