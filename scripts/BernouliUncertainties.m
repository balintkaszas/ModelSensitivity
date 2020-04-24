addPath;

ndofs = [4, 8, 16, 32, 64, 24, 48];
timeSpan = [0, pi];
Ntimes = 100;
dT = pi/Ntimes;
time2 = 0:dT:pi;
%Output = zeros(length(ndofs),length(time2));
%scaling = zeros(length(ndofs), 1);
epsilon = 1e-4;

for i=1:length(ndofs)
    n = ndofs(i);
    ndof = n;  % Even number for the beam
    ndof_spv = 2*ndof;

    ne = ndof/2;
    [M,C,K,fphi,fnl]=L_Bernoulli_Beam_Model(ne);

    % ----- First order form -----
    % \dot{x} = Ax + Fnl(x) + \epsilon*Fphi
    A = [zeros(ndof),eye(ndof);-M\K,-M\C];
    Fnl = [zeros(ndof,1);-M\fnl];
    Fphi = [zeros(ndof,1);M\fphi];
    var = symvar(Fnl);
    Fnn = matlabFunction(Fnl, 'vars', var);
    maxPlace = zeros(ndof_spv, 1).';
    lDerivative = @(t, x, e, eov) bernoullibeam(t, x, e, eov, A, Fnn, Fphi, ne);
    tic;
    [t, uncertestimate] = computeAlongTrajectoryVariationGeneralFullmodel(lDerivative, maxPlace, timeSpan, dT, 'finiteDifference');
    timeitTook = toc;
    disp(n);
    filename = sprintf('Bernoulli_uncerts_n%d.mat', n);
    save(filename, 'ndof', 't', 'epsilon', 'uncertestimate', 'timeitTook');
end



