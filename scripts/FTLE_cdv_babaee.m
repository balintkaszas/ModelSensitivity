resolution = [251,251];
domain = [0.7, 1.3;-1.,-0.4];


initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 6);
initGrid(:,1) = initialPosition(:,1);
initGrid(:,4) = initialPosition(:,2);


%charney de vore, with the 'modified' EOM
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;
timeInterval = [0,70];
stepSize = 0.1;
figure(1);
lDerivative = @(t,x,e,eov) d_charneydevore_babaee(t, x,p,eov, e);


[lmax,~] = computeCGInvariants(lDerivative, initGrid, timeInterval, 'eoV');
ftle = 0.5*log(lmax)/timeInterval(2);

x1i = linspace(0.7, 1.3,251);
x4i = linspace(-1., -0.4,251);

colormap(jet);
FF = repmat(reshape(ftle, [251, 251]),1,1,3);
%FF=smooth3(FF,'gaussian',3);
surf(x1i,x4i,reshape(ftle, [251, 251]));shading interp; axis equal;axis tight;colorbar; title('FTLE'); xlabel('z_1'); ylabel('z_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
[f_min f_max] = caxis;

