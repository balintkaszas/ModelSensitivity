domain = [-1.5, 1.5;-1.5,1.5];
resolution = [1000,1000];
initialPosition = initialize_ic_grid(resolution,domain);

%%Calculate integral (sqrt(\Lambda_s^t(x_s))) in forward time, by following
%%trajectories.
lDerivative = @(t,x,~)d_phi(t,x,0); %Duffing with dissipation 
T = 10*2*pi;
deltaT = T/1000;
ftl = uncertEstimate(lDerivative, initialPosition, resolution, [0, T], deltaT );
save('data/duffing_0_10T_forward.mat', 'ftl');
imagesc(domain(1,:),domain(2,:),ftl);
colorbar;
title('Model uncertainty (dissipative system) over $[0, 10T]$, displayed on the plane of $x_0$','Interpreter','latex')
saveas(gcf, 'pics/duffing_0_10T_forward.png');

clf;
cla;
domain = [-2, 2;-2,2];
resolution = [1000,1000];
initialPosition = initialize_ic_grid(resolution,domain);

%%Calculate integral (sqrt(\Lambda_s^t(x_s))) in forward time, by following
%%trajectories.
lDerivative = @(t,x,~)d_phi_conservative(t,x,0); %Duffing with dissipation 
T = 10*2*pi;
deltaT = T/1000;
ftl = uncertEstimate(lDerivative, initialPosition, resolution, [0, T], deltaT );
save('data/conservative_0_10T_forward.mat', 'ftl');
imagesc(domain(1,:),domain(2,:),ftl);
colorbar;
title('Model uncertainty (dissipative system) over $[0, 10T]$, displayed on the plane of $x_0$','Interpreter','latex')
saveas(gcf, 'pics/conservative_0_10T_forward.png');
