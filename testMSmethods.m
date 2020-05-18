f = @(t,x) d_abc(t,x,0, false);

resolution = [5,5]; %% coarse resolution
domain = [0, 2*pi;0,2*pi];
r  = 2;
init = Grid(3, [1,2], resolution, domain, 1e-8);

% define system:
ABC = DynSystem(f, 3, [1,1], @EOVABC);


tolerance = 1e-12;
dT = 0.1;
timeSpan = [0, 8];
disp('Regular MS starting')
tic;
MSOde45 = modelSensitivity(ABC, init, timeSpan, false, false, 0);
toc
disp('Reduced order MS (2nd order) starting')
tic;
MSOTD = modelSensitivity(ABC, init, timeSpan, false, true, 2);
toc


x1i = linspace(0, 2*pi,5);
x4i = linspace(0, 2*pi,5);
subplot(2,2,1)
title('ODE45 MS')
surf(x1i,x4i,MSOde45);shading interp; axis equal;axis tight;colorbar;  xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar

subplot(2,2,2)
title('OTD MS')
surf(x1i,x4i,MSOTD);shading interp; axis equal;axis tight;colorbar; xlabel('X_1'); ylabel('X_4')
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