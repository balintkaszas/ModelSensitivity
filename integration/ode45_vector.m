

function yf = ode45_vector(odefun,tspan,y0,useEoV)
% Reshape m-by-2 array to 8n column array
y0 = transpose(y0);
y0 = y0(:);
if useEoV
    coupledSize = 6;
else
    coupledSize = 6;
end
% Specify three timesteps for ode45's tspan. This has been reported to
% reduce memory usage.
tsteps = [tspan(1),tspan(1) + .5*diff(tspan),tspan(2)];

[~,yf] = ode45(odefun,tsteps,y0, odeset('relTol',1e-12));
yf = yf(end,:);

yf = transpose(reshape(yf,coupledSize,size(yf,2)/coupledSize));
end

