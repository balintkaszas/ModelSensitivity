%%% General case
addPath;
l = load('Amtx.mat');
Amtx = l.Amtx;
Amtx(2,1) = 0.7226;
Amtx(1,3) = 2.4248;
b = ones(3,1);
Ainvb = Amtx\b;
sigma = [0, 0, 0.8660; 0, 0.8660, 0; 0.8660, 0, 0.8660];

f0 = @(x) funzero(x, Amtx);
x0=eye(3);
x = fsolve(f0,x0)
sigma = x;


function tr = funzero(x, A)
    invSymma = inv(A.'+ A);
    tr = trace(x*x.'*invSymma) - trace(x*x.')*trace(invSymma);
end