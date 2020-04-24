function uncert = uncertSecondOrder(derivative, initialPosition, resolution, timeInterval, stepSize)
grid = initialPosition;
range = timeInterval(1):stepSize:timeInterval(2)-stepSize;
n = length(range);
%disp(n);
aggregate = zeros(size(initialPosition,1), size(range,1));
for i = 1:n
    ftle = SecondOrderComputeOverTimespan(derivative, grid, resolution, 1e-8, [range(i), timeInterval(2)]); %%Set Delta = 1e-8, the grid spacing for the aux grid.
    grid = ode45_vector(derivative, [range(i), range(i)+stepSize], grid, false);
    temp = ftle; %%integrate the square root
    %disp([iv, iv2, temp]);
	aggregate(:,i) = temp; 
end
%disp(size(range));
%disp(size(aggregate));

summ = max(aggregate,[], 2);
uncert = reshape(summ,fliplr(resolution));
end
