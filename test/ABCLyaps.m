resolution = [251,251];
%domain = [-1.5,-0.9;-1, -0.4];
domain = [0,2*pi;0, 2*pi];
range = diff(domain(2,:));
delta = range./resolution(1);
% Y-Z plane
initialPosition = initialize_ic_grid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 3);
initGrid(:,2) = initialPosition(:,1); 
initGrid(:,3) = initialPosition(:,2);

lDerivative = @(t,x, e, eov) d_abc(t, x, e,eov);

timeInterval = [0, 8];
figure(1);
title('EoV')
[ftl, ~] = computeCGInvariants(lDerivative, initGrid, timeInterval, 'eoV');

imagesc(domain(:,1), domain(:,2), log(reshape(ftl, [251, 251]))/(2*timeInterval(2)));

figure(2);
[ftl2, ~] = computeCGInvariants(lDerivative, initGrid, timeInterval, 'finiteDifference');
title('Finite Difference')
imagesc(domain(:,1), domain(:,2), log(reshape(ftl2, [251, 251]))/(2*timeInterval(2)));
