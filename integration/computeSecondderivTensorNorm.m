%%secondderiv

function norms = computeSecondderivTensorNorm(finalPosition,finalPositionEdge, finalPositionOriginal, initialPosition,resolution)
%compute_cgStrain   Compute Cauchy-Green strain
%
%   finalPosition gives positions grouped in rows. In 2 dimensions, each
%   row has 8 columns with the format:
%      (x+delta, y) (x-delta, y) (x, y+delta) (x, y-delta)
%   
%   cgStrain gives the Cauchy-Green strain tensor. In 2 dimensions it has
%   the format: c_11 c_12 c_22
% On main grid:
% finalPosition = [x1 y1; x2 y2; ... xN yN]
%
% On auxiliary grid:
% finalPosition = [x1+deltaX y1 x1-deltaX y1 x1 y1+deltaY x1 y1-deltaY;
%                  x2+deltaX y2 x2-deltaX y2 x2 y2+deltaY x2 y2-deltaY;
%                                          ...
%                  xN+deltaX yN xN-deltaX yN xN yN+deltaY xN yN-deltaY]
% finalX = x1 + deltaX; x1 - deltaX; x1;
% finalY = y1 + deltaY; y1 - deltaY; y1;
% finalPositionEdge = x1+deltaX y1+deltaY x1+deltaX y1-deltaY x1-deltaX
% y1+deltaY x1-deltaX y1-deltaY;
%finalPositionOriginal = x,y 
        norms = nan(length(initialPosition), 1);
        finalX = finalPosition(:,1:2:7);
        finalY = finalPosition(:,2:2:8);
        finalEdgeX = finalPositionEdge(:,1:2:7);
        finalEdgeY = finalPositionEdge(:,2:2:8);
        
        deltaX = diff(initialPosition(1,[3,1]))/2.; 
        deltaY = diff(initialPosition(1,[8,6]))/2.;
        
        T111 = (finalX(:,1) + finalX(:,2) - 2*finalPositionOriginal(:,1))/(deltaX.^2);
        T112 = (finalEdgeX(:,1) + finalEdgeX(:,4) - finalEdgeX(:,3) - finalEdgeX(:,2))/(4*deltaX*deltaY);
        T121 = T112;
        T122 = (finalX(:,3) + finalX(:,4) - 2*finalPositionOriginal(:,1))/(deltaY.^2);
        
        T211 = (finalY(:,1) + finalY(:,2) - 2*finalPositionOriginal(:,2))/(deltaX.^2);
        T212 = (finalEdgeY(:,1) + finalEdgeY(:,4) - finalEdgeY(:,3) - finalEdgeY(:,2))/(4*deltaX*deltaY);
        T221 = T212;
        T222 = (finalY(:,3) + finalY(:,4) - 2*finalPositionOriginal(:,2))/(deltaY.^2);

        norms = sqrt(T111.^2 + T112.^2 + T121.^2 + T122.^2 + T211.^2 + T221.^2 + T212.^2 + T222.^2);
end
