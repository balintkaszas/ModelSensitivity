addPath;
%% Define the system
dyref = @(t,x) d_duffing(t,x); %% specify derivative for the duffing system

timeSpan = [0, 2*pi]; %% time span is one period of the driving
duff = DynSystem(dyref, 2, [1,1]); % define the DynSystem object, with both Delta being 1

%% Create the computational grid
resolution = [100, 100]; 
domain = [-1.5, 1.5; -1.5, 1.5];
init = Grid(2, [1,2], resolution, domain, 1e-3);

%% Create parallel pool with 2 workers
pool = parpool('local', 2);

%% set a tolerance of 1e-7 for integration and use finite differences
ms = modelSensitivity(duff, init, timeSpan, true, 1e-7, 'finitedifference');

%set up grid for displaying the MS field
xs = linspace(domain(1,1), domain(1,2), resolution(1));
ys = linspace(domain(2,1), domain(2,2), resolution(2));

colormap(jet)
FF = repmat(ms,1,1,3);
FF=smooth3(FF,'gaussian',3);
ms1 = FF(:,:,1);
surf(xs,ys,ms);shading interp; axis equal;axis tight;colorbar; xlabel('x'); ylabel('y')
view([0 0 1]); axis equal; axis tight; shading interp;camlight;
set(gca,'colorscale', 'log');
set(gca,'ZScale', 'log');

pool.delete();