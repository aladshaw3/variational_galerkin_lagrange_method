classdef NeumannBC < BoundaryCondition    
    properties(Access = protected)
        duo     % Expected derivative value of u at boundary 
    end
    
    methods
        function obj = NeumannBC(duo)
            obj = obj@BoundaryCondition();
            obj.duo = duo;
        end
        
        function res = eval_dLdc(obj, bs, c, r, lambda)
            arguments
                obj (1,1) {mustBeA(obj,'NeumannBC')}
                bs (1,1) {mustBeA(bs,'BasisSet')}
                c (:,1) {mustBeNumeric}
                r (:,:) {mustBeNumeric}
                lambda (1,1) {mustBeNumeric}
            end  
            res = bs.dphi(r).*lambda;
        end

        function res = eval_dLdl(obj, bs, c, r, lambda)
            arguments
                obj (1,1) {mustBeA(obj,'NeumannBC')}
                bs (1,1) {mustBeA(bs,'BasisSet')}
                c (:,1) {mustBeNumeric}
                r (:,:) {mustBeNumeric}
                lambda (1,1) {mustBeNumeric}
            end
            bs.set_coeff(c);
            res = bs.du(r) - obj.duo;
        end
    end
end

