domain = [-1.5, 1.5;-1.5,1.5];
resolution = [40,40];
initialPosition = initialize_ic_grid(resolution,domain);

%%Calculate integral (sqrt(\Lambda_s^t(x_s))) in forward time, by following
%%trajectories.
lDerivative = @(t,x,~)d_phi(t,x); %Duffing with dissipation 
T = 3*2*pi;
deltaT = T/50;
ftl = uncertEstimateBackwards(lDerivative, initialPosition, resolution, [0, T], deltaT );
