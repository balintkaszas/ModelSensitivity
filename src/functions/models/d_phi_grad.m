% Forced-damped Duffing oscillator used with aperiodic forcing
function derivative_ = d_phi_grad(tau,x)

    dux = 0;
    duy = 1;
    dvx = 1-3.*x(1).^2;
    dvy = -0.15;
    derivative_ = [dux, duy; dvx, dvy];
end
