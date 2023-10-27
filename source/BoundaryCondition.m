classdef BoundaryCondition < handle
    properties (Access = protected)
    end
    
    methods
        function obj = BoundaryCondition()
        end
    end

    methods (Abstract)
        eval_dLdc(obj, bs, c, r, lambda) % Evaluate the residual of dL/dcj
        eval_dLdl(obj, bs, c, r, lambda) % Evaluate the residual of dL/dlk
    end
end

