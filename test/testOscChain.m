%%TEsting oscillator chain, against duffing oscillator.


resolution = [50,50];
domain = [-1.5, 1.5; -1.5, 1.5];
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 3);
initGrid(:,1) = initialPosition(:,1);
initGrid(:,3) = initialPosition(:,2);
n = 1;
m = 1;
k = -1;
c = 0.15;
g = 0.5;
T = 1;
alpha = 0.3;

deriv = @(t,x,e, eov) OscillatorChain(t, x, n, m, k, c, g, T, alpha, e, eov);
deriv2 = @(t,x,e, eov) d_phi(t, x, e, eov);
%deriv(1, initialPosition(1, :).')
%deriv2(1, initialPosition)
[eigmax, ~] = computeCGInvariants(deriv, initialPosition, [0, 8], 'finiteDifference');
[eigmax2, ~] = computeCGInvariants(deriv2, initialPosition, [0, 8], 'finiteDifference');
figure(1);
title('Osc Chain');
imagesc(reshape(log(eigmax), fliplr(resolution)));
colorbar;

figure(2);
title('Duffing');
imagesc(reshape(log(eigmax2), fliplr(resolution)));
colorbar;


figure(3);
title('error');
imagesc(reshape(eigmax2-eigmax, fliplr(resolution)));
colorbar;
