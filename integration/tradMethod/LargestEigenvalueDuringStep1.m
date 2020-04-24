function ft = LargestEigenvalueDuringStep1(derivative, initialPosition, timeInterval)
timelength = timeInterval(2) - timeInterval(1);
T = timelength/10;
temporaryPos = initialPosition;
finalDists = 0.;
for i = 1:10
    interval = [timeInterval(1) + (i-1)*T, timeInterval(1) + i*T];
    lderiv = @(t, x) derivative(t, x, 0, false);

    [~,yfref] = ode45(lderiv, interval, temporaryPos); 
    [~,yf] = ode45(lderiv, interval, perturbIc(temporaryPos,1e-14));
    yfref = yfref(end, :);
    yf = yf(end, :);
    
    temporaryPos = yfref;
    diff =  yf -yfref;
    dist = sqrt(sum(diff.^2,2));
    disp(dist)
    finalDists = finalDists + log(dist/1e-14);
end

ft = finalDists/timelength;
end