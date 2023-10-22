classdef PhysicsResidual < handle
    properties (Access = protected)
        lb          % lower bound of domain in which the physics apply
        ub          % upper bound of domain in which the physics apply
    end
    
    methods
        function obj = PhysicsResidual(lb, ub)
            arguments
                lb (1,:) {mustBeNumeric} = 0
                ub (1,:) {mustBeNumeric} = 1
            end
            obj.lb = lb;
            obj.ub = ub;
        end
    end

    methods (Abstract)
        eval(obj, bs_rhs, bs_lhs) % Evaluate the residual based on basis sets
    end
end

