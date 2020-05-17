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


times = linspace(0.5, 160, 20);
ms = zeros(size(times));
msanal = zeros(size(times));

for i=1:length(times)
    disp(i);
    prop = expm(A*times(i));
    
    MSSQ = modelSensitivityGlobal(deriv, maxPlace.', [1,1], [times(1), times(1) + i*DT], DT/2, 'finiteDifference', false, [1,1]);
    ms(i) = MSSQ;
    msanal(i) = compms([times(1), times(1) + i*DT], DT/2, A);
end
hold on;
plot(times, ms, 'o-');

plot(times, msanal, 'o-');
save('ms_zero_oscchain.mat', 'times', 'ms', 'msanal');

function ms = compms(timespan, dt,  A)
    times = timespan(1):dt:timespan(2)-dt;
    n = length(times);
    CGMaxToSum = zeros(size(times,1)); %% predefine the matrix containting CG maximal eigenvalues, for each subinterval and each gridpoint
    CGTraceToSum = zeros(size(times,1)); %% same for the CG trace
    for i=1:n
        T = timespan(2) - times(i);
        df = expm(A*T);
        CGMaxToSum(i) = max(svd(df));
        CGTraceToSum(i) = sum(svd(df).^2);
    end
    summCGMax = trapz(times, CGMaxToSum);
    summCGTrace = trapz(times, CGTraceToSum);
    ms = summCGMax.^2 + summCGTrace;
end
