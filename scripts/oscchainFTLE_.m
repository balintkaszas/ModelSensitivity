%%TEsting oscillator chain, against duffing oscillator.
%Duffing's parameter setup, considerable driving
addPath;
n = 16;
m = 1;
k = 1;
c = 1;
g = 0.5;
T = 1;
alpha = 1e-5;
maxPlace = zeros(32, 1);

deriv = @(t,x) OscillatorChain(t, x, n, m, k, c, g, T, alpha, 0, false);

temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));
C = c * temp;
K = k * temp;
M = m*eye(n,n);
A = [zeros(n),eye(n);-M\K,-M\C];


times = linspace(0.5*pi, 150*pi, 30);
ftle = zeros(size(times));
ftleanal = zeros(size(times));

for i=1:length(times)
    disp(i);
    prop = expm(A*times(i));
    
    [eigmax, ~] = computeCGInvariants(deriv, maxPlace.', [0, times(i)], 'finiteDifference', false);
    ftle(i) = log(eigmax)/times(i);
    ftleanal(i) = log(max(svd(prop)))/times(i);
end
hold on;
plot(times, ftleanal, 'o-');

plot(times, 0.5*ftle, 'o-');
ftle = 0.5*ftle;
save('ftle_zero_oscchain.mat', 'times', 'ftle', 'ftleanal');


