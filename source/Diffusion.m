classdef Diffusion < PhysicsResidual
    properties (Access = protected)
        D   % Diffusion coefficient
    end
    
    methods
        function obj = Diffusion(lb, ub, val)
            arguments
                lb (1,:) {mustBeNumeric} = 0
                ub (1,:) {mustBeNumeric} = 1
                val (1,1) {mustBePositive} = 1
            end
            obj = obj@PhysicsResidual(lb, ub);
            obj.D = val;
        end
        
        % bs_rhs -> j
        % bs_lhs -> i
        % return dR/dcj of diffusion
        %       = -D * sum(i,  c_i * phi_j * d2phi_i)
        function res = eval(obj, bs_rhs, bs_lhs, c_i)
            arguments
                obj (1,1) {mustBeA(obj,'Diffusion')}
                bs_rhs (1,1) {mustBeA(bs_rhs,'BasisSet')}
                bs_lhs (1,1) {mustBeA(bs_lhs,'BasisSet')}
                c_i (:,1) {mustBeNumeric}
            end
            res = -obj.D .* (integral_x_d2phi(bs_rhs, bs_lhs, obj.lb, obj.ub) * c_i);
        end
    end
end

