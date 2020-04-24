function uncert = uncertEstimateBackwards(derivative, initialPosition, resolution, timeInterval, stepSize)
range = (timeInterval(2)-stepSize:-stepSize:timeInterval(1));
n = length(range);
aggregate = zeros(size(initialPosition,1), size(range,1));
for i = 1:n
    ftle = computeOverTimespan(derivative, initialPosition, resolution, 1e-8, [timeInterval(2), range(i) ]); %%Set Delta = 1e-8, the grid spacing for the aux grid.
    temp = sqrt(ftle); %%integrate the square root
	aggregate(:,i) = temp; 
    disp(i);
end
summ = trapz(range, aggregate, 2);
uncert = reshape(summ,fliplr(resolution));
end
