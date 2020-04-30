function dy = OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov)
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
%alpha = 0.3;
f0 = ones(n,1); % loading shape
f = @(t,T)alpha*f0*cos(t/T); % loading function
A = [zeros(ndof),eye(ndof);-M\K,-M\C]; %% this is correct, just transforms to first order form.
%disp(A);
% x(1:ndof) : x
% x(ndof+1:end) : dot(x).
S = @(x)g*SP_nonlinearity(x);

Fnl = @(x) [zeros(ndof,1);-M\S(x(1:ndof))];  
Fphi = @(t, T) [zeros(ndof,1);M\f(t,T)]; %this is the time dep. forcing, only acts on velocities

dy = oscChain(t, x, Fphi, Fnl, A, T, e);
end

function dy = oscChain(t, x, Fphi, Fnl, A, T, epsilon)
    %the actual differential equation in first order form
    dy = A*x + Fnl(x) + Fphi(t, T) + epsilon;
end

function f = SP_nonlinearity(x)
%cubic nonlinear coupling. Padded with zeros for the edges
xp = [0; x; 0];
s = (xp(2:end) - xp(1:end-1)).^3;

f = s(1:end-1)-s(2:end);
end

    