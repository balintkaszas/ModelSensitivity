

function yf = ode45_vector(odefun,tspan,y0)
% Reshape m-by-2 array to 8n column array
yflat = y0(:);

% Specify three timesteps for ode45's tspan. This has been reported to
% reduce memory usage.
tsteps = [tspan(1),tspan(1) + .5*diff(tspan),tspan(2)];

[~,yf] = ode15s(odefun,tsteps,yflat, odeset('relTol',1e-12, 'absTol', 1e-13));
yf = yf(end,:);

yf = transpose(reshape(yf,size(y0)));
end

