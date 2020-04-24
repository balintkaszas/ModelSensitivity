% Forced-damped Duffing oscillator used with aperiodic forcing
function dPhi = d_phi_conservative(tau,x)
idx1 = 1:2:size(x,1);
idx2 = 2:2:size(x,1);
dPhi = nan(size(x));

dPhi(idx1) = x(idx2);
dPhi(idx2) = x(idx1) - x(idx1).^3 + 0.08*cos(tau); 

