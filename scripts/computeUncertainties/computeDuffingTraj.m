domain = [-1.5, 1.5;-1.5,1.5];
resolution = [50,50];
initialPosition = initialize_ic_grid(resolution,domain);

%%Calculate integral (sqrt(\Lambda_s^t(x_s))) in forward time, by following
%%trajectories.
ic1 = [0.8604, -0.527];%, 3791.];
 ic2 = [0.9835, -0.47];%, 128.]
T = 3*2*pi;
deltaT = T/100;
%ip1 = [0.2147, 0.3408];
timeInterval = [0, T];
stepSize = deltaT;
lDerivative = @(t,y) d_phi(t,y,0, false);
[t1, uncertAtTime1] = computeAlongTrajectory(lDerivative, ic1, timeInterval, stepSize );
[t2, uncertAtTime2] = computeAlongTrajectory(lDerivative, ic2, timeInterval, stepSize );
save('uncertAtTime1.mat', 'uncertAtTime1');
save('t1,mat', 't1');
save('uncertAtTime2.mat', 'uncertAtTime2');
save('t2.mat', 't2');



% 

%ftl = uncertEstimate(lDerivative, initialPoint, [1,1], [0, T], stepSize );
%imagesc(domain(:,1), domain(:,2), ftl);
%colorbar;
