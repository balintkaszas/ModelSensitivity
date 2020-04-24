% Forced-damped Duffing oscillator used with aperiodic forcing

function dPhi = d_linear(tau,x, e, useEoV, alpha)
if useEoV
    idx1 = 1:6:size(x,1)-5;
    idx2 = 2:6:size(x,1)-4;
else
    idx1 = 1:2:numel(x)-1;
    idx2 = 2:2:numel(x);
end
dPhi = nan(size(x));

dPhi(idx1) = alpha*x(idx1);
dPhi(idx2) = alpha*x(idx2); 
if useEoV
    % Define terms of the equation of variation
    idx3 = 3:6:size(x,1)-3;
    idx4 = 4:6:size(x,1)-2;
    idx5 = 5:6:size(x,1)-1;
    idx6 = 6:6:size(x,1);
    
    dux = alpha;
    duy = 0;
    dvx = 0;
    dvy = alpha;
    
    % Perform matrix multiplication manually
    dPhi(idx3) = dux.*x(idx3) + duy.*x(idx5);
    dPhi(idx4) = dux.*x(idx4) + duy.*x(idx6);
    dPhi(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
    dPhi(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
end

