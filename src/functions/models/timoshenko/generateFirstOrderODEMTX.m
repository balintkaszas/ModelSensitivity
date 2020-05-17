function A = generateFirstOrderODEMTX(system)
%from the system class that contains symbolic functions describing the
%second order system, it returns the right hand side of the first order
%system and the gradient of the system.

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

% compute the gradient symbolically 
Fgrad = jacobian(Fnl, X);

Fnn = matlabFunction(Fnl, 'vars', X);
FgradNum = matlabFunction(Fgrad, 'vars', X);

A = double(A);
Asparse = sparse(A);
odeFunction = @(t, xvar, e, eov) odfunc(t, xvar, e, eov, A, Fnn);
gradOdeFunction = @(t, xvar, e, eov) gradfunc(t, xvar, e, eov, Asparse, FgradNum);

end

function odf = odfunc(t, xvar, e, eov, A, Fnn)
    unpacked = num2cell(xvar);
    odf = A*xvar + Fnn(unpacked{:}); % Fnn expects separate variables as argument, instead of a vector..
end


function odf = gradfunc(t, xvar, e, eov, A, FgradNum ) %returns the jacobian matrix A+\nabla Fnl(X)
    unpacked = num2cell(xvar);
    odf = A ;%+ FgradNum(unpacked{:}); % FgradNum expects separate variables as argument, instead of a vector..
end
