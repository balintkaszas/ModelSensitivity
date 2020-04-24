% Forced-damped Duffing oscillator used with aperiodic forcing
function derivative_ = d_phi_conservative(tau,x, epsilon, useEoV)
if useEoV
    idx1 = 1:6:size(x,1)-5;
    idx2 = 2:6:size(x,1)-4;
else
    idx1 = 1:2:numel(x)-1;
    idx2 = 2:2:numel(x);
end

derivative_ = nan(size(x));

derivative_(idx1) = x(idx2);
derivative_(idx2) = x(idx1) - x(idx1).^3 + 0.08*cos(tau) + epsilon;
if useEoV
    % Define terms of the equation of variation
    idx3 = 3:6:size(x,1)-3;
    idx4 = 4:6:size(x,1)-2;
    idx5 = 5:6:size(x,1)-1;
    idx6 = 6:6:size(x,1);
    
    dux = 0;
    duy = 1;%1-3.*x(idx1).^2;
    dvx = 1-3.*x(idx1).^2;
    dvy = 0.;
    %%FIXME
    % Perform matrix multiplication manually
    derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
    derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
    derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
    derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
end
