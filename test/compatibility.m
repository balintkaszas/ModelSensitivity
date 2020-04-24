%load forward and backward uncertainties
b = load('../data/areBAndFCompatible/back_3_cons.mat');
back = b.ftl;
f = load('../data/areBAndFCompatible/forw_3_cons.mat');
forw = f.ftl2;
%% regular grid: 
domain = [-1.5, 1.5;-1.5,1.5];
resolution = [1000,1000];
initialPosition = initialize_ic_grid(resolution,domain,2);
%derivatives and timeframe:
lDerivative = @(t,x,~)d_phi(t,x, 0, false); %Duffing with dissipation 
T = 3*2*pi;
deltaT = T/50;
%% U(x_0) = U(F_t^t0(x_t))
%RHS:
advectedGrid = ode45_vector(lDerivative, [T, 0], initialPosition, false);
%%Interpolate back to a regular grid, but now at x_0! 
%https://ch.mathworks.com/help/matlab/math/interpolation-using-a-specific-delaunay-triangulation.html
DT = delaunayTriangulation(advectedGrid);
vi = nearestNeighbor(DT,initialPosition);
Vq = back(vi);
figure(1);
hold on;
imagesc(domain(1,:),domain(2,:), reshape(abs(Vq), fliplr(resolution)));
colorbar;
set(gca,'YDir','normal');
set(gca,'colorscale','log');

%plot(advectedGrid(:,1), advectedGrid(:,2), '.');
hold off;
figure(2);
imagesc(domain(1,:),domain(2,:),forw);
set(gca,'YDir','normal');
set(gca,'colorscale','log');
colorbar;

