classdef Monomial1D < BasisSet

    methods (Access = public)
        %% Constructor for Monomial1D
        function obj = Monomial1D(N)
            arguments
                N (1,1) {mustBePositive} = 1
            end
            obj = obj@BasisSet(N);
            obj.Nqp = 2.*obj.Nqp-1;
        end

        %% Basis function for Monomial1D
        %
        %       Calling this function will evaluate the basis functions
        %       for a given set of positions. Since this is a 1D function,
        %       the shape of position vector r is Rx1 and we calculate 
        %       each phi as...
        %           phi = r^(i-1)   for i = [1, N]
        %
        %   @param r Rx1 spatial positions
        %
        %   @return phi NxR matrix of N basis functions in R positions
        function phi = phi(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'Monomial1D')}
                r (:,1) {mustBeNumeric}
            end
            phi = (r.^((1:obj.N)-1))'; % return NxR vec of each func evaluated at r
        end

        %% 1st Derivative of Basis function for Monomial1D
        %
        %       Calling this function will evaluate the 1st derivatives
        %       for a given set of positions. Since this is a 1D function,
        %       the shape of position vector r is Rx1 and we calculate 
        %       each phi as...
        %           dphi = (i-1)*r^(i-2)   for i = [1, N]
        %
        %   @param r Rx1 spatial positions
        %
        %   @return dphi NxR matrix of N 1st derivatives in R positions
        function dphi = dphi(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'Monomial1D')}
                r (:,1) {mustBeNumeric}
            end
            dphi = (r.^((1:obj.N)-2))'.*((1:obj.N)'-1); % return NxR vec of each func evaluated at r 
            dphi(isnan(dphi))=0; % Correct for error near r=0
        end

        %% 2nd Derivative of Basis function for Monomial1D
        %
        %       Calling this function will evaluate the 2nd derivatives
        %       for a given set of positions. Since this is a 1D function,
        %       the shape of position vector r is Rx1 and we calculate 
        %       each phi as...
        %           d2phi = (i-1)*(i-2)*r^(i-3)   for i = [1, N]
        %
        %   @param r Rx1 spatial positions
        %
        %   @return d2phi NxR matrix of N 2nd derivatives in R positions
        function d2phi = d2phi(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'Monomial1D')}
                r (:,1) {mustBeNumeric}
            end
            d2phi = (r.^((1:obj.N)-3))'.*((1:obj.N)'-1).*((1:obj.N)'-2); % return NxR vec of each func evaluated at r
            d2phi(isnan(d2phi))=0; % Correct for error near r=0
        end
    end
end