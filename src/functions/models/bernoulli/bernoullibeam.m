function df = bernoullibeam(tau,x, epsilon, useEoV, A, Fnl, Fphi, n)
    eps = 0.01;
    %size(Fnl)
    df = A*x + Fnl(x(2*n-1,1)) + Fphi;
end
