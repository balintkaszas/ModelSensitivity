
%% When the difference between the bound and the uncertainty is positive, the inequality fails

epsilon = .001;

time = linspace(0,2*pi, 51);
deriv = @(t,y) d_phi(t,y,0, false);
[timm, uncerr] = computeAlongTrajectorySq(deriv, transpose(mnPlace), [0., 2*pi], 2*pi/50);
sensit = uncerr;

f = @(t,x) d_phi(t,x,0, false);
g = @(t,x) epsilon;
dt = 1e-5*2*pi;
time = linspace(0,2*pi, 51);
time2 = 0:dt:2*pi;

[time2,yfref] = ode45(@(t,y) d_phi(t,y,0, false),time2,mnPlace); 
aggregate = zeros(size(yfref,1),100);
hold on;

for i = 1:100
    yf = sde_euler(f,g,time2,mnPlace);
    diff = yf-yfref;
    aggregate(:,i) = diff(:,1).^2 + diff(:,2).^2;
    plot(time2, aggregate(:,i), '-.');
    
    %plot(t, sqrt(diff(:,1).^2+diff(:,2).^2), '.');
end
plot(time2, mean(aggregate, 2), '-', 'linewidth', 3, 'color', 'black');
plot(time, 4*squeeze(sensit)*epsilon.^2, '-o', 'linewidth', 3,'color', 'red');
xlabel('Time, $t$', 'interpreter', 'latex');
ylabel('$||x(t,\varepsilon)-x(t,0)||^2$','interpreter', 'latex');
legend(gca,  '$|x(t) - x_0(t)|^2$', '$4\varepsilon^2\int_{t_0}^{t_1}\Lambda_s^{t_1}(x_s)ds$', 'interpreter', 'latex')
title('$\varepsilon = 0.001$', 'interpreter', 'latex');
