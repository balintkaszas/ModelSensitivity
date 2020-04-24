function [x,y] = detectRidge(scalarField, resolution, domain)
%% Ridge detection following T. Lindeberg (1998)


    %computing the gradient:

    %grid spacing in either direction
    differenceX = diff(domain(1,:))/resolution(1);
    differenceY = diff(domain(2,:))/resolution(2);


    %generating the uniform grid
    initialPosition = initialize_ic_grid(resolution, domain, 2);
    coords = reshape(initialPosition, [resolution(1), resolution(2), 2]);
    xi = coords(:,:,1);
    yi = coords(:,:,2);


    %calculate the gradient
    [Lx,Ly] = gradient(scalarField, differenceX, differenceY);

    % calculate Hessian
    [Lxx, Lxy] = gradient(Lx, differenceX, differenceY);
    [Lyx, Lyy] = gradient(Ly, differenceX, differenceY);

    detShapeOp = sqrt((Lxx-Lyy).^2 + 4.*Lxy.^2);

    %calculate rotation angles: 
    cosb = sqrt(0.5 * (1 + (Lxx - Lyy)./detShapeOp));
    sinb = sign(Lxy).*sqrt(0.5 * (1 - (Lxx - Lyy)./detShapeOp));

    %calculate directional derivatives
    Lp = sinb.*Lx - cosb.*Ly;
    [Lpx,Lpy] = gradient(Lp, differenceX, differenceY);
    [Lqx,Lqy] = gradient(Lp, differenceX, differenceY); %%d/dx(dL/dp)

    Lpp = sinb.*Lpx - cosb.*Lpy;
    Lqq = cosb.*Lqx + sinb.*Lqy;
    ridgeMask = contour(xi, yi, Lp, [0,0]);
    xCandidate = ridgeMask(1,:);
    yCandidate = ridgeMask(2,:);
    xi = xi';
    yi = yi';
    Lpp = Lpp';
    %interpolants
    Lppinterp = griddedInterpolant(xi,yi,Lpp, 'linear');
    Lqqinterp = griddedInterpolant(xi,yi,Lqq, 'linear');

    ridgeMask = (Lppinterp(xCandidate, yCandidate) < 0) & (abs(Lppinterp(xCandidate, yCandidate)) >= abs(Lqqinterp(xCandidate, yCandidate))); %% ridge mask satisfying two conditions
    x = xCandidate(ridgeMask);
    y = yCandidate(ridgeMask);
    
end