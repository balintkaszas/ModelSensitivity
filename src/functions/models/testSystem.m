% Forced-damped Duffing oscillator used with aperiodic forcing
function derivative_ = testSystem(tau,x, epsilon, useEoV)

% disp(x(idx1));
% disp(x(idx2));
r = sqrt(x(1).^2 + x(2).^2);
derivative_ = nan(size(x));

derivative_(1) = x(1) - x(1).^3;% - (1./r - r).*(x(1).^2 - x(2).^2).*x(2)).*cos(2*r) + (x(2) - 0.5*x(2)*(3*x(1).^2 + x(2).^2) - (1./r -r ).*2*x(1).*x(2).^2).*sin(2*r);
derivative_(2) = -x(2);%: + x(2).^3;% + (1./r - r).*(x(1).^2 - x(2).^2).*x(1)).*cos(2*r) + (x(1) - 0.5*x(1)*(3*x(2).^2 + x(1).^2) - (1./r -r ).*2*x(2).*x(1).^2).*sin(2*r);
if(useEoV)
    dux = 1 - 3*x(1).^2;
    duy = 0;
    dvx = 0;
    dvy = -1;% + 3*x(2).^2;
    derivative_(3) = dux.*x(3) + duy.*x(5);
    derivative_(4) = dux.*x(4) + duy.*x(6);
    derivative_(5) = dvx.*x(3) + dvy.*x(5);
    derivative_(6) = dvx.*x(4) + dvy.*x(6);
end
