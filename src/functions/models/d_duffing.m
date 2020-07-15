% Forced-damped Duffing oscillator used with periodic forcing
function derivative_ = d_duffing(tau,x)

derivative_ = nan(size(x));

derivative_(1) = x(2);
derivative_(2) = x(1) - x(1).^3 - 0.15*x(2) + 0.3*cos(tau);
end
