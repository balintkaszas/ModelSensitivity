% Forced-damped Duffing oscillator used with aperiodic forcing
function derivative_ = d_pendulum(tau, x, epsilon, useEoV, tramp, cmax)

idx1 = 1:2:numel(x)-1;
idx2 = 2:2:numel(x);

derivative_ = nan(size(x));

derivative_(1) = x(2);
derivative_(2) = -0.1*x(2) - (1./9.)*sin(x(1)) + Cramp(tau, tramp, cmax)*cos(x(1))*cos(tau) + epsilon;

end



function C = Cramp(tau, tramp, Cmax)

C0 = 0.1;
m = (Cmax-C0)/(tramp*2*pi);
C = C0;
if tau >= 0 && tau <= tramp*2*pi
    C=C0 + tau*m;
else
    C = C0 + m*tramp*2*pi;
end
end

