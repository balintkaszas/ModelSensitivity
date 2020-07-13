

testGridClass()


function passes = testGridClass()
asd = Grid(3, [1,2],[100,100], [0,2*pi; 0,2*pi], 1e-5);

%compare against 'manual method'
resolution = [100,100];
domain = [0,2*pi;0, 2*pi];
% X-Z plane
initialPosition = initializeGrid(resolution, domain, 2);
initGrid = zeros(length(initialPosition), 3);
initGrid(:,1) = initialPosition(:,1); 
initGrid(:,2) = initialPosition(:,2);

passes = all(all(initGrid == asd.points));
end