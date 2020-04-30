function dy = grad_d_oscChain(t, x, n, m, k, c, g, T, alpha)
%Returns the RHS of an oscillator chain model.
%n: dimension
%m: uniform mass
%k: linear coupling constant
%c: damping coeff
%g: nonlinearity coeff
%T: forcing period -> T=1 means period of 1.
%alpha: forcing stuff

% m = 1;
% k = 1;
% c = 1;
% g = 0.5;

temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));
temp(1,1) = 1;
C = c * temp;
K = k * temp;
M = m*eye(n,n);


ndof = n;
%disp(A);
% x(1:ndof) : x
% x(ndof+1:end) : dot(x).
S = @(x)g*SP_nonlinearity_derv(x);

gradmtx = @(x) -M\S(x(1:ndof));
dy = [zeros(ndof),eye(ndof);-M\K + gradmtx(x),-M\C]; %% this is correct, just transforms to first order form.

end



function Df = SP_nonlinearity_derv(x)
n = length(x);
X = [0; x; 0];
s = 3*(X(2:end) - X(1:end-1)).^2;

D = [-s(2:end-1); 0]; 
E = s(1:end-1)+s(2:end);
F = [0; -s(2:end-1)]; 

Df = spdiags([D E F], -1:1, n,n);
end
    