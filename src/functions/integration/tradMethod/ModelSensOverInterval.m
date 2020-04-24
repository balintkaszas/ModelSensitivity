function uncert = ModelSensOverInterval(derivative, initialPosition, resolution, timeInterval, stepSize)
grid = initialPosition;
range = timeInterval(1):stepSize:timeInterval(2) - stepSize;
n = length(range);
%disp(n);
aggregate = zeros(size(initialPosition,1), size(range,1));
for i = 1:n
    ftle = LargestEigenvalueDuringStep1(derivative, grid, [range(i), timeInterval(2)]); %%implicit renormalization in each step
    lderiv = @(t, x) derivative(t, x, 0, false);

    grid = ode45_vector(lderiv, [range(i), range(i) + stepSize], grid, false);
    %temp = sqrt(ftle); %%integrate the square root
    %disp([iv, iv2, temp]);
    T = timeInterval(2)-range(i); %%t-s
   	aggregate(:,i) = ftle; 
end
summ = trapz(range, aggregate, 2);
uncert = reshape(summ,fliplr(resolution));
end

