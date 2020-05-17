
p.C = 0.1;
p.z1Star = 0.95;
p.z4Star = -0.76095;
p.beta = 1.25;
p.gamma = 0.2;
p.b = 0.5;

knull = zeros(6,6);
f = @(t,x) d_abc(t,x,0, false);

resolution = [250,250];
domain = [0, 2*pi;0,2*pi];
r  = 2;
init = Grid(3, [1,2], resolution, domain, 1e-8);
tolerance = 1e-12;
dT = 0.1;
timemax = 8;
MM = computeOTD_ODE45(f, @EOVABC, init.points, [0, timemax],r, dT);
%save('FTLETimoshenko_GOOOD_big_8.mat', 'ftle', 'timeSpan', 'resolution', 'domain', 'vars');
ftle = reshape(log(MM)/(2*timemax), resolution);
figure(2)
colormap(jet)
%FF = repmat(ftle_otd,1,1,3);
%FF=smooth3(FF,'gaussian',3);
x1i = linspace(0, 2*pi,251);
x4i = linspace(0, 2*pi,251);
surf(x1i,x4i,ftle);shading interp; axis equal;axis tight;colorbar; title(['FTLE OTD(r=', num2str(r),')']); xlabel('X_1'); ylabel('X_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
colorbar
% for i = 1:length(times)
%     timemax = times(i);
%     dT = times(i)/20;
%     disp(i);
%     mm = computeOTD_ODE45(f, @EOVABC, initialPosition, [0, timemax],r, dT);
%     ftle(i) = log(mm)/(2*timemax);
% end
% 
% plot(times, ftle, '-');


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