%%% Dissipative case
a = load('data/stable.mat');
b = load('data/unstable.mat');
u = load('uncertanity/data/good_uncert_dissipative_10T__asd.mat');
stableManifold = a.finalPositionMainGridBack;
unstableManifold = b.finalPositionMainGridD;
uncert = u.ftl;
mnPlace = [0.3889, -0.521];
hold on;
domain = [-1.5,1.5;-1.5, 1.5];
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
set(gca,'colorscale','log');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\dot{x}$', 'Interpreter', 'latex');
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])
%plot(mPlace(1),mPlace(2),'.','color', 'Red', 'MarkerSize',25);
%plot(mnPlace(1),mnPlace(2),'.','color', 'Green', 'MarkerSize',25);
%saveas(gcf, 'pics/dissipative_markPoints.png');
cla;
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
set(gca,'colorscale','log');
leg = legend('Stable manifold','Unstable manifold');
set(leg, 'Interpreter','latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\dot{x}$', 'Interpreter', 'latex');
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])
saveas(gcf, 'pics/dissipative_.png');
cla;

%%% Conservative case
a = load('data/stableConservative.mat');
b = load('data/unstableConservative.mat');
u = load('uncertanity/data/good_uncert_conservative_10T__asd.mat');
stableManifold = a.finalPositionMainGridBack;
unstableManifold = b.finalPositionMainGridD;
uncert = u.ftl;
domain = [-2,2;-1.5, 1.5];
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
set(gca,'colorscale','log');
leg = legend('Stable manifold','Unstable manifold');
set(leg, 'Interpreter','latex');
xlabel('x', 'Interpreter', 'latex');
ylabel('$\dot{x}$', 'Interpreter', 'latex');

xlim([-2, 2])
ylim([-1.5, 1.5])
saveas(gcf, 'pics/conservative_.png');
%cla;
%[1.103, -0.1186];
%81.99

%[0.9349, -0.8904]
%2064
% 
% 
% %%% Positive lyapunov case
% a = load('data/stablePoslyap.mat');
% b = load('data/unstablePoslyap.mat');
% u = load('uncertanity/data/good_uncert_poslyap_10T__asd.mat');
% stableManifold = a.finalPositionMainGridBack;
% unstableManifold = b.finalPositionMainGridD;
% uncert = u.ftl;
% domain = [-2,2;-1.5, 1.5];
% imagesc(domain(1,:),domain(2,:),uncert);
% colorbar;
% plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
% plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
% set(gca,'colorscale','log');
% leg = legend('Stable manifold','Unstable manifold');
% set(leg, 'Interpreter','latex');
% xlabel('x', 'Interpreter', 'latex');
% ylabel('$\dot{x}$', 'Interpreter', 'latex');
% 
% xlim([-2, 2])
% ylim([-1.5, 1.5])
% saveas(gcf, 'pics/poslyap_.png');
