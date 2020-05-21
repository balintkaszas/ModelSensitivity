%initialize_ic_grid Initialize initial conditions on cartesian grid
function position = initialize_ic_grid(resolution,domain)

dimGrid = size(domain, 1);
position=zeros(prod(resolution), dimGrid);

    if dimGrid == 2
        xVector = linspace(domain(1,1),domain(1,2),resolution(1));
        yVector = linspace(domain(2,1),domain(2,2),resolution(2));
        [positionX,positionY] = meshgrid(xVector,yVector);
        position(:,1) = positionX(:);
        position(:,2) = positionY(:);
    elseif dimGrid == 3
        xVector = linspace(domain(1,1),domain(1,2),resolution(1));
        yVector = linspace(domain(2,1),domain(2,2),resolution(2));
        zVector = linspace(domain(3,1),domain(3,2),resolution(3));

        [positionX,positionY, positionZ] = meshgrid(xVector,yVector, zVector);
        position(:,1) = positionX(:);
        position(:,2) = positionY(:);
        position(:,3) = positionZ(:);
    else 
        disp("Only 2D or 3D grids possible")
    end

end

