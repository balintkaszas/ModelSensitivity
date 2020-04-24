


n = 2;
temp = spdiags(ones(n,1)*[-1 2 -1], [-1 0 1], sparse(n,n));
%temp(1,1) = 1;
C = 1 * temp;
K = 1 * temp;
M = 1*eye(n,n);
A = [zeros(n),eye(n);-M\K,-M\C]; %% this is correct, just transforms to first order form.
full(K)

%sd = svd(A);
%plot(asd.^2)