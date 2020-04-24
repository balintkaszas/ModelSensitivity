resolution = [250,250];
domain = [-pi, pi; -pi, pi];
initialPosition = initialize_ic_grid(resolution, domain, 2);


lDerivative = @(t,x,e,eov) d_pendulum(t, x,eov, e, 3, 0.33);

pool = parpool('local', 16);

[helo,~] = computeCGInvariants(lDerivative, initialPosition, [0, 15*2*pi], 'finiteDifference', false);
imagesc(reshape(log(helo), resolution));
saveas(gcf,'Barchart.png');
save('TIPPING3_33.mat', 'helo');
pool.delete();