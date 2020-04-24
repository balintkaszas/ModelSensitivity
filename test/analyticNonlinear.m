function result = analyticNonlinear(pos, T, dt)
%% Validate the numerics on a nonlinear model. 
%% The analytic expression for \Lambda^t_s(xs) is maxx(x0, t, tmax) 

range = (0:dt:T-dt);
aggregate = zeros(size(range,1));
n = length(range);
for i = 1:n
    time = range(i);
    ft = maxx(pos(1), time, T);
    temp = sqrt(ft);
    aggregate(i) = temp;
end
result = trapz(range, aggregate);
%uncert = reshape(summ,fliplr(resolution)); 
%result = max(max((uncert - ftl)./uncert)); %% Error is the maximal relative deviation from the analytic expression (over the grid)
end

function maxeigen = maxx(x0, t, tmax)
lambda = exp(tmax-t);
x = x0*exp(t);
x11 = lambda^2;
x12 = 2*x*(lambda^3-1)/3;
x22 = (4/9)*x.^2*(lambda^3-1)^2/lambda^2 + 1/lambda^2;
maxeigen = (x11 + x22 + sqrt((x11-x22).^2 + 4*x12.*x12))/2; % explicit formula for 2by2 eigenvalue
end

