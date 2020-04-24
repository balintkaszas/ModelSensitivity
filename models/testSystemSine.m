function derivative_ = testSystemSine(tau,x, epsilon, useEoV)

derivative_ = nan(size(x));

derivative_(1) = x(1) - sin(x(2)) - (x(1) - sin(x(2))).^3 - x(2).*cos(x(2));
derivative_(2) = -x(2);
if(useEoV)
    dux = 1 - 3*(x(1) - sin(x(2))).^2;
    duy = -cos(x(2)) + 3*(x(1) - sin(x(2))).^2.*cos(x(2)) + x(2).* sin(x(2)) - cos(x(2));
    dvx = 0;
    dvy = -1;
    derivative_(3) = dux.*x(3) + duy.*x(5);
    derivative_(4) = dux.*x(4) + duy.*x(6);
    derivative_(5) = dvx.*x(3) + dvy.*x(5);
    derivative_(6) = dvx.*x(4) + dvy.*x(6);
end
