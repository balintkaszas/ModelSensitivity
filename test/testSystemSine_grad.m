function dy = testSystemSine_grad(tau,x)

    dux = 1 - 3*(x(1) - sin(x(2))).^2;
    duy = -cos(x(2)) + 3*(x(1) - sin(x(2))).^2.*cos(x(2)) + x(2).* sin(x(2)) - cos(x(2));
    dvx = 0;
    dvy = -1;

    dy = [dux, duy;dvx, dvy];
end

