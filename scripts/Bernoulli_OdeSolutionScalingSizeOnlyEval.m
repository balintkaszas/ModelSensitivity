cla;clf;
addPath;
%size scaling:
%%No ODE solution, only function evaluation. My guess is that ODE solution takes O(n) function evaluations, and the function evaluations are O(n) themselves.
sizes = 2*floor(linspace(2, 3000, 15));
S = zeros(3, length(sizes));
for i = 1:length(sizes)
    ndof = sizes(i);  % Even number for the beam
    A =  randi([0, 5], [ndof,ndof]);
    b = randi([0, 5], [ndof,1]);
    lDerivative = @(t,x) withMatmul(t, x, A, b);  
    lDerivative2 = @(t,x) withoutMatmul(t, x, b);   

	disp(sizes(i))
%	tic
    x0 = randi([0, 5], [ndof,1]);
	f = @() lDerivative(1, x0);
	g = @() lDerivative2(1, x0);

    S(1, i) = sizes(i);
	S(2, i) = timeit(f);
    S(3, i) = timeit(g);
    
end
save('SizeScalingEvals.mat', 'S');

function dx = withMatmul(t, x, A, b)
    dx = A*x + b;
end

function dx = withoutMatmul(t, x, b)
    dx = b;
end
