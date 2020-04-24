
function ftleT = computeOverTimespan(derivative, initialPosition, resolution, auxGridRelDelta, timespan)
           %% Eigenvectors from auxiliary grid
            initialPositionM = reshape(initialPosition,[fliplr(resolution),2]);
            deltaX = mean(diff(initialPositionM(1,:,1)));
            deltaY = mean(diff(initialPositionM(:,1,2)));

            gridSpace = deltaX;
            if(size(initialPosition, 1) ==1 )
                gridSpace = 1.;
            end
            auxiliaryGridAbsoluteDelta = gridSpace*auxGridRelDelta;
            auxiliaryPosition = auxiliary_position(initialPosition,auxiliaryGridAbsoluteDelta);

            % Transform auxiliaryPosition into a two column array
            auxiliaryPositionX = auxiliaryPosition(:,1:2:end-1);
            auxiliaryPositionY = auxiliaryPosition(:,2:2:end);
            auxiliaryPositionC = [auxiliaryPositionX(:),auxiliaryPositionY(:)];
            
            finalPositionAuxGrid = ode45_vector(@(t,y)derivative(t,y),timespan,auxiliaryPositionC, false);

            % Transform finalPosition into an eight column array
            finalPositionAuxGridX = finalPositionAuxGrid(:,1);
            finalPositionAuxGridY = finalPositionAuxGrid(:,2);
            nPoints = prod(double(resolution));
            finalPositionAuxGridX = reshape(finalPositionAuxGridX,nPoints,4);
            finalPositionAuxGridY = reshape(finalPositionAuxGridY,nPoints,4);
            finalPositionAuxGrid = nan(nPoints,8);
            finalPositionAuxGrid(:,1:2:7) = finalPositionAuxGridX;
            finalPositionAuxGrid(:,2:2:8) = finalPositionAuxGridY;


            cgStrainAuxGrid = compute_cgStrain(finalPositionAuxGrid,auxiliaryPosition,resolution);
            [cgStrainV,cgStrainD] = arrayfun(@(x11,x12,x22)eig_array(x11,x12,x22),cgStrainAuxGrid(:,1),cgStrainAuxGrid(:,2),cgStrainAuxGrid(:,3),'UniformOutput',false);
                cgStrainD = cell2mat(cgStrainD);
            %cgStrainD = %negative_to_nan(cgStrainD);
            cgEigenvalue2 = (cgStrainD(:,2));

            ftleT = cgEigenvalue2;%ftle(cgEigenvalue2,diff(timespan));
end
