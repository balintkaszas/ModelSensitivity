% derivative Double gyre velocity field
%
% SYNTAX
% derivative_ = derivative(t,position,useEoV,epsilon,amplitude,omega)
%
% INPUT ARGUMENTS
% t: time
% position: [x1;y1;x2;y2;...;xn;yn]
% useEov: logical that controls use of the equation of variation
% epsilon,amplitude,omega: double gyre parameters
%
% REFERENCE
% DOI:10.1016/j.physd.2005.10.007

function derivative_ = derivative2(t,x,e,useEoV,epsilon,amplitude,omega)

validateattributes(t,{'double'},{'scalar'})
validateattributes(x,{'double'},{'column'})

validateattributes(x,{'double'},{'column'})
validateattributes(useEoV,{'logical'},{'scalar'})
validateattributes(epsilon,{'double'},{'scalar'})
validateattributes(amplitude,{'double'},{'scalar'})
validateattributes(omega,{'double'},{'scalar'})

a = epsilon*sin(omega*t);
b = 1 - 2*epsilon*sin(omega*t);
forcing = a*x(1).^2 + b*x(1);

derivative_ = nan(size(x));

derivative_(1) = -pi*amplitude*sin(pi*forcing).*cos(pi*x(2));
derivative_(2) = pi*amplitude*cos(pi*forcing).*sin(pi*x(2)).*(2*a*x(2) + b);
if(useEoV)
    dux = -pi^2*amplitude*cos(pi*forcing).*cos(pi*x(2)).*(2*a*x(1) + b);
    duy = pi^2*amplitude*sin(pi*forcing).*sin(pi*x(2));
    dvx = -pi^2*amplitude*sin(pi*forcing).*sin(pi*x(2)).*(2*a*x(1) + b) + 2*a*pi*amplitude*cos(pi*forcing).*sin(pi*x(2));
    dvy = pi^2*amplitude*cos(pi*forcing).*cos(pi*x(2)).*(2*a*x(1) + b);

    % Perform matrix multiplication manually
    derivative_(3) = dux.*x(3) + duy.*x(5);
    derivative_(3) = dux.*x(4) + duy.*x(6);
    derivative_(5) = dvx.*x(3) + dvy.*x(5);
    derivative_(6) = dvx.*x(4) + dvy.*x(6);
end