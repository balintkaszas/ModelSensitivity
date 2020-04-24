classdef Grid
    %GRID Represents the computational grid
    %   Detailed explanation goes here
    
    properties
        dimension 
        domain
        resolution
        variables
        points
        difference
    end
    
    methods
        function obj = Grid(dimension, variables, resolution, domain, difference)
            %GRID Construct an instance of this class
            %   Detailed explanation goes here
            obj.dimension = dimension;
            obj.domain = domain;
            obj.variables = variables;
            obj.resolution = resolution;
            obj.difference = difference;
            nonZeroPoints = initialize_ic_grid(resolution, domain, 2);
            obj.points = zeros(size(nonZeroPoints,1), dimension);
            for i = 1:length(variables) % fill in the places of the nonzero variables
                obj.points(:, variables(i)) = nonZeroPoints(:,i);
            end

        end
    end
end

