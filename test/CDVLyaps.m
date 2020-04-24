%%Charney de vore, tikme dependent FTLE from specific IC
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;
timeInterval = [0,70];
lDerivative = @(t,x,e,eov) d_charneydevore_babaee(t, x,p,eov, e);



for i=1:length(times)
    time = times(i);
    helo = computeOverTimespanVariationalNonvectorized(lDerivative, IC, [0,time]);
    lyaps(i,1) = 0.5*log(helo(1,1))/time;
    lyaps(i,2) = 0.5*log(helo(1,2))/time;
    lyaps(i,3) = 0.5*log(helo(1,3))/time;
end
hold on;
plot(times, lyaps(:,1), '.');
plot(times, lyaps(:,2), '.');
plot(times, lyaps(:,3), '.');
xlabel('$T$', 'interpreter', 'latex');
ylabel('FTLE$(\Lambda_0^T)$', 'interpreter', 'latex');
ylim([-0.2, 0.65]);
% figure(2);
% lDerivative = @(t,x) d_abc(t, x, false);
% helo2 = computeOverTimespanFiniteDiffGeneral(lDerivative, initGrid,delta,1e-8, [0,8]);
% ftle2 = reshape(log(helo2)/8, [251,251]);
% imagesc(domain(1,:), domain(2,:), ftle2);
% title('Finite Difference');
% set(gca, 'Ydir', 'normal');
% colorbar;
% 
% figure(3);
% imagesc(domain(1,:), domain(2,:), abs(ftle2-ftle)./ftle);
% set(gca, 'Ydir', 'normal');
% colorbar;
% title('Rel. Error');

%%CORRECT 
% [t,helo] = LargestEigenvalueDuringStep1(lDerivative, [1.14, 0., 0, -0.91,0.,0.], [0,100]);
% plot(t, helo)



%ModelSensOverInterval(lDerivative, initGrid, resolution, timeInterval, stepSize)
%helo = ode45_vector(lDerivative, [0, 100], initGrid, false);


