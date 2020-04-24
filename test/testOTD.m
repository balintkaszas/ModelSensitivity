
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;

knull = zeros(6,6);
f = @(t,x,e , eov) d_abc(t,x,e, eov);

resolution = [250,250];
domain = [0, 2*pi;0,2*pi];

initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 3);
initGrid(:,1) = initialPosition(:,1);
initGrid(:,3) = initialPosition(:,2);
%tic
%[eigmax, ~] = computeCGInvariants(f, initGrid, [0, 8], 'eoV');
%toc
%tic
%eigotd = computeOTD(f, initGrid, [0, 8], 2, 8./100);
%toc
 Ag= load('ABC_Flow_ftle.mat');
 FTLE_exact = Ag.FTLE_exact;
 FTLE_OTD = Ag.FTLE_OTD;
 
kmax = size(FTLE_exact,2);
ftle_exact = reshape(FTLE_exact(:,kmax),251,251);
ftle_otd = reshape(FTLE_OTD(1,:),251,251);
x1i = linspace(0, 2*pi,251);
x4i = linspace(0, 2*pi,251);
figure(1);
colormap(jet)
%FF = repmat(ftle_exact,1,1,3);
%FF=smooth3(FF,'gaussian',3);
surf(x1i,x4i,ftle_exact);shading interp; axis equal;axis tight;colorbar; title('FTLE(Exact)'); xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
[f_min f_max] = caxis;
colorbar;
r=2;
figure(2)
colormap(jet)
%FF = repmat(ftle_otd,1,1,3);
%FF=smooth3(FF,'gaussian',3);
surf(x1i,x4i,ftle_otd);shading interp; axis equal;axis tight;colorbar; title(['FTLE OTD(r=', num2str(r),')']); xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar


figure(3);
colormap(jet)
surf(x1i,x4i,-(ftle_otd-ftle_exact));shading interp; axis equal;axis tight;colorbar; title(['FTLE OTD(r=', num2str(r),')']); xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar