function yf = ode45Ensemble(fun, tspan,y0, isParallel)
%% Wrapper to advect a grid of points, with a non-vectorized derivative

nRows = size(y0,1);
yf = zeros(size(y0)); 

% Add a single intermediate step in timespan, to not save all
% intermedaite values. This saves memory.
tsteps = [tspan(1), (tspan(1)+tspan(2))/2, tspan(2)];

%% Separate case for a single trajectory
if nRows == 1
    y = y0(1,:);
    y = y(:);
    [~,solution] = ode15s(fun, tsteps, y, odeset('relTol',1e-12));  
    solution = solution(end,:);
    yf(1,:) = transpose(solution);    
end
%% Loop over gridpoints
if nRows > 1
    if isParallel == true %use parallel for if possible
        parfor i=1:nRows
            y = y0(i,:);   %get current grid point
            y = y(:);  %flatten it
            [~,solution] = ode15s(fun, tsteps, y, odeset('relTol',1e-12));  
            solution = solution(end,:);
            yf(i,:) = transpose(solution);
        end
    else
        for i=1:nRows %use regular for otherwise
            y = y0(i,:);
            y = y(:);
            [~,solution] = ode15s(fun,tsteps,y, odeset('relTol',1e-12));  
            solution = solution(end,:);
            yf(i,:) = transpose(solution);
        end
    end
end
end
