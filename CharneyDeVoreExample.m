%% Define the system
%CDV params:
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;


dyref = @(t,x) d_charneyDeVore(t,x,p, false);
dygrad = @(t,x) d_charneyDeVore_grad(t,x,p);
cdv = DynSystem(dyref, 6, [1,1], dygrad); % 6 dimensional system

%% Create the computational grid
resolution = [250, 250];
domain1 = [-3, 3; -3, 3];
init = Grid(6, [2,4], resolution, domain1, 1e-8); %Z_2 - z_4
timeSpan = [0, 15];

%% Create parallel pool with 2 workers
pool = parpool('local', 2);

%% set a tolerance of 1e-7 for integration and use finite differences
ms = modelSensitivity(cdv, init, timeSpan, true, 1e-7, 'finitedifference');

%set up grid for displaying the MS field
xs = linspace(domain(1,1), domain(1,2), resolution(1));
ys = linspace(domain(2,1), domain(2,2), resolution(2));

colormap(jet)
FF = repmat(reshape(ms, [250, 250]),1,1,3);
FF=smooth3(FF,'gaussian',3);
ms1 = FF(:,:,1);
surf(xs,ys,ms1);shading interp; axis equal;axis tight;colorbar; xlabel('z_2'); ylabel('z_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight;
set(gca,'colorscale', 'log');
set(gca,'ZScale', 'log');

pool.delete();