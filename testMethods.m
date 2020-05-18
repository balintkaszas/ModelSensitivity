f = @(t,x) d_abc(t,x,0, false);

resolution = [50,50]; %% coarse resolution
domain = [0, 2*pi;0,2*pi];
r  = 2;
init = Grid(3, [1,2], resolution, domain, 1e-8);

% define system:
ABC = DynSystem(f, 3, [0,0], @EOVABC);


tolerance = 1e-12;
dT = 0.1;
timeSpan = [0, 8];
disp('Regular FTLE starting')
tic;
[FTLEODE45, TRACEODE45] = FTLE(ABC, init, timeSpan, false, 1e-8, false, 0);
toc
disp('Reduced order FTLE (2nd order) starting')
tic;
[FTLEOTD, TRACEOTD] = FTLE(ABC, init, timeSpan, false, 1e-8, true, 2);
toc
disp('TRUE FTLE starting')
tic;
[FTLEoriginal, TRACEoriginal] = computeCGInvariants(f, init.points, timeSpan, 'finiteDifference', false, 1e-8);
toc
x1i = linspace(0, 2*pi,50);
x4i = linspace(0, 2*pi,50);

subplot(3,2,1)
title('ODE45 FTLE')
surf(x1i,x4i,FTLEODE45);shading interp; axis equal;axis tight;colorbar;  xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar

subplot(3,2,2)
title('ODE45 TRACE')
surf(x1i,x4i,TRACEODE45);shading interp; axis equal;axis tight;colorbar; xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar



subplot(3,2,3)
title('OTD FTLE')
surf(x1i,x4i,FTLEOTD);shading interp; axis equal;axis tight;colorbar;xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar


subplot(3,2,4)
title('OTD TRACE')
surf(x1i,x4i,TRACEOTD);shading interp; axis equal;axis tight;colorbar;  xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar


subplot(3,2,5)
title('Control FTLE')
surf(x1i,x4i,reshape(FTLEoriginal, resolution));shading interp; axis equal;axis tight;colorbar; xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar


subplot(3,2,6)
title('Control TRACE')
surf(x1i,x4i,reshape(TRACEoriginal, resolution));shading interp; axis equal;axis tight;colorbar;  xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar


function eov = EOVABC(t, x)
A = sqrt(3);
B = sqrt(2);
C = 1;
    J11 = 0;
    J12 = -C*sin(x(2));
    J13 = A * cos(x(3));
    
    J21 = B*cos(x(1));
    J22 = 0;
    J23 = -A*sin(x(3));
    
    J31 = -B*sin(x(1));
    J32 = C*cos(x(2));
    J33 = 0;
    eov = [J11, J12, J13; J21, J22, J23; J31, J32, J33];
end