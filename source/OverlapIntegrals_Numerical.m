classdef OverlapIntegrals_Numerical
    methods (Static)

        %% Function to compute the overlap integrals between variable 
        %       j and variable i between lb and ub bounds. The object 
        %       will pick the appropriate number of qp points based on 
        %       each individual variables' requirements. 
        %
        %       overlap = integral{ var_j'(r) * fun_i(r) * dV }
        %
        %   @param var_j Right-hand side variable object
        %   @param fun_i Left-hand side function object 
        %   @param lb 1xd vector of lower bounds (where d is dimension)
        %   @param ub 1xd vector of upper bounds (where d is dimension)
        function overlap = phi_x_fun(var_j, fun_i, lb, ub)
            arguments
                var_j (1,1) {mustBeA(var_j,'BasisSet')}
                fun_i (1,1) {mustBeA(fun_i,'function_handle')}
                lb (1,:) {mustBeNumeric}
                ub (1,:) {mustBeNumeric}
            end
            if (length(ub) ~= length(lb))
                ME = MException('OverlapIntegrals_Numerical:invalidSize',...
                    '\nInput vectors must be of same size/dimension, d_lb=%g, d_ub=%g',string(length(lb)),string(length(ub)));
                throw(ME)
            end

            Nqp = max(var_j.get_Nqp())*100;
            pts = zeros(Nqp,length(lb));
            pts(1,:) = lb;
            pts(end,:) = ub;
            for i=1:length(lb)
                pts(:,i) = lb(i):(ub(i)-lb(i))/(Nqp-1):ub(i);
            end
            dr = (ub-lb)./(Nqp-1);
            dV = prod(dr);
            overlap = (0.5*fun_i(pts(1:end-1,:))*var_j.phi(pts(1:end-1,:))'*dV + ...
                        0.5*fun_i(pts(2:end,:))*var_j.phi(pts(2:end,:))'*dV)';
        end

    end
end