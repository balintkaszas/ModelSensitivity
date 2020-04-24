function ic2 = perturbIc(ic, delta)

%randomvector = -1 + 2.*rand(size(ic));
randomvector = randn(size(ic));
randomvector = bsxfun(@rdivide,randomvector,sqrt(sum(randomvector.^2,2)));
ic2 = ic+delta.*randomvector;
end
