cla;clf;
addPath;
%size scaling:
sizes = 2*floor(linspace(2, 200, 8));
S = zeros(2, length(sizes));
for i = 1:length(sizes)
    ndof = sizes(i);  % Even number for the beam
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
    
    f = @(t,x,e,eov) bernoullibeam(t, x, 0, false, A, Fnn, Fphi, ne);
    disp(sizes(i))
    maxPlace = zeros(ndof_spv, 1);
    tic;
    [~, ~] = computeCGInvariants(f, maxPlace, [0, 2*pi], 'finiteDifference');
    S(1, i) = sizes(i);
	S(2, i) = toc;
end
save('SizeScalingWithEuler.mat', 'S');




