cla;clf;
addPath;




ndof = 4;  % Even number for the beam
ndof_spv = 2*ndof;

ne = ndof/2;
[M,C,K,fphi,fnl]=L_Bernoulli_Beam_Model(ne);

% ----- First order form -----
% \dot{x} = Ax + Fnl(x) + \epsilon*Fphi
A = [zeros(ndof),eye(ndof);-M\K,-M\C];
Fnl = [zeros(ndof,1);-M\fnl];
Fphi = [zeros(ndof,1);M\fphi];
%[t, sol] = ode45(@(t,x) bernoullibeam(t, x, 0, false, A, Fnl, Fphi), [0,10], zeros(20,1));
resolution = [10,10];
Fnn = matlabFunction(Fnl, 'vars', 'u_3');

domain = [-1.5, 1.5;-1.5, 1.5]*1e-3;

%initialPosition = initialize_ic_grid(resolution, domain, 2);
%initgrid = zeros(length(initialPosition), ndof_spv);
%initgrid(:, ndof) = initialPosition(:, 1);
%initgrid(:, 2*ndof) = initialPosition(:, 2);

%%Time scaling:
	times = linspace(0.1, 100, 15);
	lDerivative = @(t,x) bernoullibeam(t, x, 0, false, A, Fnn, Fphi);
T = zeros(2, length(times));
for i = 1:length(times)
	disp(times(i))
	tic
	[~,sol] = ode45(lDerivative, [0,times(i)], zeros(ndof_spv,1));
	T(1, i) = times(i);
	T(2, i) = toc;
end
save('Timescaling.mat', 'T', 'ndof');




