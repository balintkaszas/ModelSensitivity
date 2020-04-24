domain = [-1.5, 1.5;-1.5,1.5];
resolution = [100,100];
timespan=[0,2*pi];
initialPosition = initialize_ic_grid(resolution,domain,2);

figure(1);
derivative = @(t,x)d_phi(t,x,0,false);
dt= 0.05*2*pi;
ftl = uncertEstimate(derivative, initialPosition, resolution, timespan, dt);
ftl = reshape(ftl, [100,100]);
imagesc(domain(1,:), domain(2,:), ftl);
title('Finite Difference');
set(gca, 'colorscale', 'log');
set(gca, 'Ydir', 'normal');
colorbar;

figure(2);
derivative = @(t,x,e,eov)d_phi(t,x,e,eov);
dt= 0.05*2*pi;
ftlV = uncertEstimateVariational(derivative, initialPosition, [100,100], timespan, dt);
ftlV = reshape(ftlV, [100,100]);
imagesc(domain(1,:), domain(2,:), ftlV);
title('Variational Equation');

set(gca, 'colorscale', 'log');
set(gca, 'Ydir', 'normal');
colorbar;

figure(3);
derivative = @(t,x)d_phi(t,x,0,false);
dt= 0.005*2*pi;
ftll = ModelSensOverInterval(derivative, initialPosition, resolution, timespan, dt);
ftll = reshape(ftll, [100,100]);
imagesc(domain(1,:), domain(2,:), ftll);
title('Only largest LE');
set(gca, 'colorscale', 'log');
set(gca, 'Ydir', 'normal');
colorbar;


figure(4);
imagesc(domain(1,:), domain(2,:), abs(ftl-ftll)./ftl);
title('Relative error: LLE/FiniteDiff');

set(gca, 'colorscale', 'log');
set(gca, 'Ydir', 'normal');
colorbar;


figure(5);
imagesc(domain(1,:), domain(2,:), abs(ftlV-ftl)./ftl);
title('Relative error: EOV/FiniteDiff');
set(gca, 'colorscale', 'log');
set(gca, 'Ydir', 'normal');
colorbar;



