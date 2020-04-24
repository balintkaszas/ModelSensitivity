%%% Dissipative case
a = load('data/stable.mat');
b = load('data/unstable.mat');
u = load('uncertanity/data/areBAndFCompatible/forw_3.mat');
stableManifold = a.finalPositionMainGridBack;
unstableManifold = b.finalPositionMainGridD;
uncert = u.ftl2;
%mnPlace = [0.3889, -0.521];
hold on;
domain = [-1.5,1.5;-1.5, 1.5];
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
set(gca,'colorscale','log');
set(gca,'YDir','normal');
%xlabel('$x$', 'Interpreter', 'latex');
%ylabel('$\dot{x}$', 'Interpreter', 'latex');
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
saveas(gcf, 'pics/oct2_dissipative_0_3T.png');
savefig('pics/oc2_dissipative_0_3T');
hold off;

cla;
clf;
domain = [-1.5,1.5;-1.5, 1.5];
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
set(gca,'colorscale','log');
set(gca,'YDir','normal');
%xlabel('$x$', 'Interpreter', 'latex');
%ylabel('$\dot{x}$', 'Interpreter', 'latex');
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
saveas(gcf, 'pics/oct2_dissipative_0_3T_nomanifold.png');
savefig('pics/oc2_dissipative_0_3T_nomanifold');
clf;
cla;
%%% Conservative case
a = load('data/stableConservative.mat');
b = load('data/unstableConservative.mat');
u = load('uncertanity/data/areBAndFCompatible/forw_3_cons.mat');
stableManifold = a.finalPositionMainGridBack;
unstableManifold = b.finalPositionMainGridD;
uncert = u.ftl2;
domain = [-1.5,1.5;-1.5, 1.5];
hold on;
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
plot(stableManifold(:,1), stableManifold(:,2),'.','MarkerSize',2, 'color','black');
plot(unstableManifold(:,1), unstableManifold(:,2),'.','MarkerSize',2,'color','red');
set(gca,'colorscale','log');
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
%xlabel('x', 'Interpreter', 'latex');
%ylabel('$\dot{x}$', 'Interpreter', 'latex');

%xlim([-2, 2])
%ylim([-2, 2])
saveas(gcf, 'pics/oct2_conservative_0_3T.png');
savefig('pics/oc2_conservative_0_3T');%cla;
hold off;
clf;
cla;
imagesc(domain(1,:),domain(2,:),uncert);
colorbar;
set(gca,'colorscale','log');
xlabel('x', 'Interpreter', 'latex');
ylabel('$\dot{x}$', 'Interpreter', 'latex');
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
%ylim([-2, 2])
saveas(gcf, 'pics/oct2_conservative_0_3T_nomanifold.png');
savefig('pics/oc2_conservative_0_3T_nomanifold');
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
% c
% xlim([-2, 2])
% ylim([-1.5, 1.5])
% saveas(gcf, 'pics/poslyap_.png');
