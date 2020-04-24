
function [cgStrainD,varargout] = negative_to_nan(cgStrainD,varargin)

negIdx = any(cgStrainD <= 0,2);

if negIdx
    warning([mfilename,':negativeEigenvalue'],'Negative eigenvalues')
end

cgStrainD(negIdx,:) = nan;

if nargout == 2
    cgStrainV = varargin{1};
    cgStrainV(negIdx,:) = nan;
    varargout{1} = cgStrainV;
end
end