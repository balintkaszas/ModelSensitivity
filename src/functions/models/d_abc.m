function dx = d_abc(tau,x,e,useEOV)
dx = nan(size(x));
A = sqrt(3);
B = sqrt(2);
C = 1;
dx(1) = A*sin(x(3)) + C* cos(x(2));
dx(2) = B*sin(x(1)) + A* cos(x(3)) + e;
dx(3) = C*sin(x(2)) + B* cos(x(1));

if useEOV
    J11 = 0;
    J12 = -C*sin(x(2));
    J13 = A * cos(x(3));
    
    J21 = B*cos(x(1));
    J22 = 0;
    J23 = -A*sin(x(3));
    
    J31 = -B*sin(x(1));
    J32 = C*cos(x(2));
    J33 = 0;
    J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
    M = reshape(x(4:end), 3,3);
    M =transpose(M);
    rm = J*M;
    rm = reshape(transpose(rm), size(x(4:end)));
    dx(4:end) = rm;
end

