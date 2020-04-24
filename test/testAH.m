%%Testing the Implementation on a transformed nonlinear saddle problem. 
domain = [-1, 1;-1,1];
resolution = [100,100];
initialPosition = initialize_ic_grid(resolution,domain,2);
A1 = initialPosition(:,1);
A2 = initialPosition(:,2);

r = sqrt(x1.^2 + x2.^2);
grid = zeros(size(initialPosition));
grid(:,1) = x1.*cos(r) - x2.*sin(r); %%transforming the ICs
grid(:,2) = x2.*cos(r) + x1.*sin(r);

a1 = reshape(grid(:,1), resolution);
a2 = reshape(grid(:,2), resolution); %%New ICs

timespan = [0, 2];
derivative = @(t,x,e,eov)testSystem(t,x,e,eov);

%% Numerical FTLE field and trace. 
[eigmaxFiniteDiff, traceFiniteDiff] = computeCGInvariants(derivative, initialPosition,  timespan, 'finiteDifference');
ftl = reshape(log(ftl)/(4), resolution);
surf(reshape(grid(:,1), resolution),reshape(grid(:,2), resolution),ftl);shading interp; axis equal;axis tight;colorbar; title('FTLE'); xlabel('z_1'); ylabel('z_4')
view([0 0 1]); axis equal; axis tight; shading interp;camlight
[f_min f_max] = caxis;


function [eigmax,tr] = analyticFTLE(A1, A2, t0, t)
    
    X1 = sign(A1).*(1 - (1-1./A1.^2)*exp(-2*(t-t0)));
    X2 = sign(A2).*(1 - (1-1./A2.^2)*exp(2*(t-t0))); %Solutions to original, separable equations
    r = sqrt(X1.^2 + X2.^2); %distance
    r0 = sqrt(A1.^2 + A2.^2);
    %the transformation is:
    %f1(X1, X2) = X1 cos(r) - X2 sin(r)
    %f2(X1, X2) = X2 cos(r) + X1 sin(r)
    
    %% Transformation of initial conditions:
    a1 = A1.*cos(r0) - A2.*sin(r0);
    a2 = A2.*cos(r0) + A1.*sin(r0);
    
    %%Differential of old ICs w.r.t new ICs.
    dA1da1 = (1 - a1.*a2./r0).*cos(r) - (a1.^2/r0).*sin(r0);
    dA1da2 = -(1 + a1.*a2./r0).*sin(r) - (a2.^2/r0).*cos(r0);
    dA2da1 = (1 - a1.*a2./r0).*sin(r) + (a1.^2/r0).*cos(r0);
    dA2da2 = (1 + a1.*a2./r0).*cos(r) - (a2.^2/r0).*sin(r0);
    
    %%Original flowmap gradient (diagonal)
    dX1dA1 = (exp(-2*(t-t0))./A1.^3).*(1 - (1 - 1./A1.^2).*exp(-2*(t-t0))).^(-3/2);  %this is the time dependent part. Cross terms are 0. 
    dX2dA2 = (exp(2*(t-t0))./A2.^3).*(1 - (1 - 1./A2.^2).*exp(2*(t-t0))).^(-3/2);  
    
    %% dX/da = dX/dA*dA/da
    dX1da1 = dX1dA1.*dA1da1;
    dX1da2 = dX2dA2.*dA1da2;
    dX2da1 = dX1dA1.*dA2da1;
    dX2da2 = dX2dA2.*dA2da2;
    drda1 = dX1da1.*X1./r + dX2da1.*X2./r;
    drda2 = dX1da2.*X1./r + dX2da2.*X2./r;
    %% transformed flowmap gradient (in new coords)
    dx1da1 = (dX1da1 - X2.*drda1).*cos(r) - (dX2da1 + X1.*drda1).*sin(r);
    dx1da2 = (dX1da2 - X2.*drda2).*cos(r) - (dX2da2 + X1.*drda2).*sin(r);
    
    
    dx2da1 = (dX2da1 + X1.*drda1).*cos(r) + (dX1da1 - X2.*drda1).*sin(r);
    dx2da2 = (dX2da2 + X1.*drda2).*cos(r) + (dX1da2 - X2.*drda2).*sin(r);
    
    %%components of the Cauchy-Green strain tensor
    C11 = dx1da1.^2 + dx2da1.^2;
    C12 = dx1da1.*dx1da2 + dx2da1.*dx2da2;
    C22 = dx1da2.^2 + dx2da2.^2;
    %% Eigenvalues of the CG tensor
    lambda1 = (C11 + C22 - sqrt((C11 - C22).^2 + 4*C12.^2))/2;
    lambda2 = (C11 + C22 + sqrt((C11 - C22).^2 + 4*C12.^2))/2;
    eigmax = lambda2;
    tr = lambda1 + lambda2;
end