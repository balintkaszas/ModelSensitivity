%%smalldfuff


addPath;

dy = @(t,x) d_phi_asd(t,x,0,false);

point = zeros(32,1);%*1e-8;
timeSpan = [0,6];
dT = 1e-7;
r = 5;

init = Grid(2,[1,2], [100, 100], [-1.5, 1.5; -1.5, 1.5], 1e-8);
[ftl, ~] = computeCGInvariants(dy, init.points, [0, 4], 'finiteDifference', false, 1e-8);
imagesc(reshape(log(ftl), [100,100]));
%time = linspace(0, derivative, initialPosition, timeSpan, method, isParallel, tolerance