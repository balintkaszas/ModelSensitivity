function [cgEigenvalueMax, cgTrace] = computeCGInvariantsWrapper(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, difference, tolerance, method)
%% Wrapper to choose correct function depending on integration method
% method = 'eov' or 'finitedifference'
% in case method = 'eov', provide arbitrary values for
%   difference
% in case method = 'finitedifference', provide an arbitrary handle for derivativeEOV
% Usage:
% [lambdaMax,Traces] =
% computeCGInvariantsWrapper(derivative,derivativeEOV, initialPosition, timeSpan, isParallel, tolerance);

    switch method
    case 'eov'
        [cgEigenvalueMax, cgTrace] = computeCGInvariantsEOV(derivative, derivativeEOV, initialPosition, timeSpan, isParallel, tolerance);
    case 'finitedifference'
        [cgEigenvalueMax, cgTrace] = computeCGInvariantsFiniteDiff(derivative, initialPosition, timeSpan, isParallel, difference, tolerance);
    otherwise
        disp('Unknown calculation method');
    end
end

