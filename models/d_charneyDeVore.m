function dx = d_charneyDeVore(tau,x, params, useEoV,e,k, omega)

dx = nan(size(x));

gamma = params.gamma;
b = params.b;
beta = params.beta;
epsilon = 16 * sqrt(2) / (5 * pi);
C = params.C;
z1Star = params.z1Star;
z4Star = params.z4Star;

dx(1) = gstar(gamma, 1, b) .* x(3) - C .* (x(1) - z1Star);
dx(2) = -(alphm(1, b) .* x(1) - betam(beta, 1, b)).* x(3) - C * x(2) - deltam(1, b) .* x(4) .* x(6);
dx(3) = (alphm(1, b) .* x(1) - betam(beta, 1, b)) .* x(2) - gamm(gamma, 1, b) .* x(1) - C*x(3) + deltam(1, b)*x(5).*x(4);
dx(4) = gstar(gamma, 2, b) .* x(6) - C*(x(4) - z4Star) + epsilon*(x(2).*x(6) - x(3).*x(5));
dx(5) = -(alphm(2, b)*x(1) - betam(beta, 2, b)).*x(6) - C * x(5) - deltam(2, b) * x(4).*x(3);
dx(6) =  (alphm(2, b)*x(1) - betam(beta, 2, b)).*x(5) - gamm(gamma, 2, b) * x(4) - C*x(6) + deltam(2,b) *x(4).*x(2);
perturbationSize = zeros(size(x));
if e > 0
    perturbationSize= e.*error(x,tau,k, omega);
end

dx = dx + perturbationSize;
if useEoV
    J11 = -C;
    J12 = 0;
    J13 = gstar(gamma, 1, b);
    J14 = 0;
    J15 = 0;
    J16 = 0;
    
    J21 = -alphm(1, b) * x(3);
    J22 = -C;
    J23 = -(alphm(1, b) * x(1) - betam(beta, 1, b));
    J24 = -deltam(1, b) * x(6);
    J25 = 0;
    J26 = -deltam(1, b) * x(4);
    
    J31 = alphm(1, b) * x(2) - gamm(gamma, 1, b);
    J32 = alphm(1, b) * x(1) - betam(beta, 1, b);
    J33 = -C;  
    J34 = deltam(1, b) * x(5);
    J35 = deltam(1, b) * x(4);
    J36 = 0;

    J41 = 0;
    J42 = epsilon * x(6);
    J43 = -epsilon * x(5);
    J44 = -C;
    J45 = -epsilon * x(3);
    J46 = gstar(gamma, 2, b) + epsilon * x(2);
    
    J51 = -alphm(2, b) * x(6);
    J52 = 0;
    J53 = -deltam(2, b) * x(4);
    J54 = -deltam(2, b) * x(3);
    J55 = -C;
    J56 = betam(beta, 2, b) - alphm(2, b) * x(1);
    
    J61 = alphm(2, b) * x(5);
    J62 = deltam(2, b) * x(4);
    J63 = 0;
    J64 = -gamm(gamma, 2, b) + deltam(2, b) * x(2);
    J65 = alphm(2, b) * x(1) - betam(beta, 2, b);
    J66 = -C;
    

    J = [J11 J12 J13 J14 J15 J16; J21 J22 J23 J24 J25 J26; J31 J32 J33 J34 J35 J36; J41 J42 J43 J44 J45 J46; J51 J52 J53 J54 J55 J56; J61 J62 J63 J64 J65 J66];
    M = reshape(x(7:end), 6,6);
    M = transpose(M);
    rm = J*M;
    rm = reshape(transpose(rm), size(x(7:end)));
    dx(7:end) = rm;
end
end

function gs = gstar(gamma, m, b)
gs = gamma * (4 * m / (4 * m^2 - 1)) * (sqrt(2) * b / pi);
end

function g = gamm(gamma, m, b)
g = gamma * (4 * m^3/(4 * m^2-1)) * (sqrt(2) * b / (pi * (b^2 + m^2)));
end

function bet = betam(beta, m, b)
bet = beta * b^2 / (b^2 + m^2);
end

function alph = alphm(m, b)
alph = (8 * sqrt(2)/pi) * (m^2 / (4 * m^2 - 1)) * ((b^2 + m^2 - 1) / (b^2 + m^2));
end

function delta = deltam(m, b)
delta = ((64 * sqrt(2))/(15 * pi)) * (b^2 - m^2 +1) / (b^2 + m^2);
end 


function errorvector = error(x, t,k, omega)
    z = 0.6; %length scale
    errorvector = zeros(6,1);
    
    for i=1:6
        ki = k(:,i)/z;
        errorvector(i) = sin(dot(ki,x))*sin(omega*t);
    end
end