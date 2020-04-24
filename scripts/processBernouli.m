%%ProcessBernouli 
%%Manual entry of running times...
runtimes = [23070, 25897, 34353, 37796, 43588];
ndofs = [4, 8, 16, 24, 32];
BB = load('Bernoulli_MC.mat');
U = load('Bernoulli_uncerts_n8.mat');
tuncert = U.t;
uncertestimate = U.uncertestimate;
Montecarlos = BB.Montecarlos;
yfref = BB.yfref;
time2 = BB.time2;
ndof = BB.ndof;
aggregate = zeros(101, 5000);
hold on;
%plot(yfref(:, 4), yfref(:, 8), '.');
rng(0, 'twister'); %%for random samples
figure(1);
hold on;
avg = averageMC(yfref, Montecarlos);
figure(2);
hold on;
for i = 1:50
    r = randi([1 5000],1);
    yf = Montecarlos(:, :, r);
    diff = yfref - yf;
    error = sum(diff.^2, 2);
    figure(1);
    plot(yf(:, ndof-1), yf(:, 2*ndof-1), '.');
    figure(2);
    
    plot(time2, error, '.');
end
figure(1);
plot(yfref(:, ndof-1), yfref(:, 2*ndof-1), '.-', 'Linewidth', 3, 'color', 'black', 'DisplayName','$\varepsilon$ = 0');
figure(2);
plot(time2, avg, '.-', 'Linewidth', 3, 'color', 'black');
plot(tuncert, uncertestimate.*1e-8*16, '.-', 'color', 'red');
%histogram(Montecarlos(101, ndof-1, :))
hold off;
figure(3);
histogram(Montecarlos(101, ndof-1, :))

function expectedVal = averageMC(yfref, Montecarlos)
    diff = yfref - Montecarlos;
    aggr = sum(diff.^2, 2);
    expectedVal = mean(aggr, 3);
end