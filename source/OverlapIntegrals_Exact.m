classdef OverlapIntegrals_Exact

    methods (Static)
        %% Function to compute the overlap integrals between Monomial1D sets 
        %       j and variable i between lb and ub bounds. The object 
        %       will pick the appropriate number of qp points based on 
        %       each individual variables' requirements. 
        %
        %       overlap = integral{ var_j.phi(r)' * var_i.phi(r) * dV }
        %
        %   @param var_j Right-hand side Monomial1D object
        %   @param var_i Left-hand side Monomial1D object 
        %   @param lb 1xd vector of lower bounds (where d is dimension)
        %   @param ub 1xd vector of upper bounds (where d is dimension)
        function overlap = phi_x_phi_Monomial1D_x_Monomial1D(var_j, var_i, lb, ub)
            arguments
                var_j (1,1) {mustBeA(var_j,'Monomial1D')}
                var_i (1,1) {mustBeA(var_i,'Monomial1D')}
                lb (1,1) {mustBeNumeric}
                ub (1,1) {mustBeNumeric}
            end

            Nj = (1:length(var_j.get_coeff()))';
            Ni = 1:length(var_i.get_coeff());
            overlap = ((ub.^(Nj+Ni-1)) - (lb.^(Nj+Ni-1)))./(Nj+Ni-1);
        end
    end
end