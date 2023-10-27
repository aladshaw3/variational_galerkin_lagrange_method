classdef DirichletBC < BoundaryCondition
    properties(Access = protected)
        uo      % Expected boundary value of the variable u
    end
    
    methods
        function obj = DirichletBC(uo)
            obj = obj@BoundaryCondition();
            obj.uo = uo;
        end
        
        function res = eval_dLdc(obj, bs, c, r, lambda)
            arguments
                obj (1,1) {mustBeA(obj,'DirichletBC')}
                bs (1,1) {mustBeA(bs,'BasisSet')}
                c (:,1) {mustBeNumeric}
                r (:,:) {mustBeNumeric}
                lambda (1,1) {mustBeNumeric}
            end  
            res = bs.phi(r).*lambda;
        end

        function res = eval_dLdl(obj, bs, c, r, lambda)
            arguments
                obj (1,1) {mustBeA(obj,'DirichletBC')}
                bs (1,1) {mustBeA(bs,'BasisSet')}
                c (:,1) {mustBeNumeric}
                r (:,:) {mustBeNumeric}
                lambda (1,1) {mustBeNumeric}
            end
            bs.set_coeff(c);
            res = bs.u(r) - obj.uo;
        end
    end
end

