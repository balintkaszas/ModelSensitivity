%% compare gradient components

A = load('A.mat');
A = A.A;

asd = @(x) TIMOSHENKODE_Grad_NL(x);

fullode = @(t,x) TIMOSHENKODE_Pert(t, x.');
time = linspace(0, 4, 100);
rng(10);
x = rand(2,1);
x0 = zeros(32,1);
x0(16) = x(1);
x0(32) = x(2);
[~, sol] = ode45(fullode, time, x0, odeset('relTol', 1e-9));
M = arrayfun(@(idx)asd(sol(idx,:)),1:100, 'UniformOutput', false); 

Ms = zeros(32,32,100);
for i = 1:100
   Ms(:,:,i) = cell2mat(M(i)); 
end
%M = cell2mat(M);
[argvalue, argmax] = max(mean(Ms, [1,2]));
disp(argvalue)
disp(argmax)
A(A==0) = nan;

MMax = Ms(:,:,argmax);
MMax(MMax ==0) = nan;

figure(1)
imagesc(A)
set(gca, 'colorscale', 'log')
colorbar;

figure(2)
imagesc(MMax)
set(gca, 'colorscale', 'log')
colorbar