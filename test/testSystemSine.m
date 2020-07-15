function derivative_ = testSystemSine(tau,x)

derivative_ = zeros(size(x));

derivative_(1) = x(1) - sin(x(2)) - (x(1) - sin(x(2))).^3 - x(2).*cos(x(2));
derivative_(2) = -x(2);
end
