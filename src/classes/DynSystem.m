classdef DynSystem
    %DYNSYSTEM Represents a dynamical system
    %   Detailed explanation goes here
    
    properties
        rhsFn
        gradRhsFn
        dimension
        deltas
    end
    
    methods
        function obj = DynSystem(rhs, dimension, deltas, varargin)
            %DYNSYSTEM construct from elementary functions
            %   Should provide the:
            %    Right-hand side of the ODE
            %    Dimension of the system (dim)
            %    Estimates for the model-uncertainty: deltas = [D1, D2] 
            %          D1 is the estimate for deterministic uncertainty,
            %          and D2 is for stochastic uncertainty.
            %  Can provide:
            %    gradient of the right hand side.
            obj.rhsFn = rhs;
            obj.dimension = dimension;
            obj.deltas = deltas;
            if(nargin > 3) % at least 3 arguments must be given 
                obj.gradRhsFn = varargin{1}; %assign the optional argument
            end
    end
        
        function dy = rhs(obj, t, x, e, eov)
            %METHOD1 Evaluate right hand side
            %   Detailed explanation goes here
            dy = obj.rhsFn(t,x,e,eov);
        end
        function dy = gradRhs(obj, t, x)
            %METHOD1 Evaluate gradient of right hand side
            %   Detailed explanation goes 
            if(obj.gradRhsFn)
                dy = obj.gradRhsFn(t, x);
            end
            
            
        end
    end
end

