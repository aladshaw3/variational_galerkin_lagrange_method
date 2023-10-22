classdef BasisSet < handle

    properties (Access = protected)
        N       % Number of basis functions
        c       % Weight coefficients for the basis functions
        Nqp     % Number of quadrature points for each basis function 
    end

    methods (Access = public)
        %% Constructor for BasisSet
        function obj = BasisSet(N)
            arguments
                N (1,1) {mustBePositive} = 1
            end
            obj.N = N;
            obj.c = zeros(N,1);
            obj.Nqp = (1:obj.N)';
        end

        %% Function to set coefficients for the BasisSet
        function [] = set_coeff(obj,c)
            arguments
                obj (1,1) {mustBeA(obj,'BasisSet')}
                c (:,1) {mustBeNumeric}
            end
            if (length(c) ~= obj.N)
                ME = MException('BasisSet:invalidSize',...
                    '\nInput vector must be of same size as "BasisSet", N=%g',string(obj.N));
                throw(ME)
            end
            obj.c = c;
        end

        %% Function to get the coefficients
        function c = get_coeff(obj)
            c = obj.c;
        end

        % Function to get the number of quadrature points for each basis
        function Nqp = get_Nqp(obj)
            Nqp = obj.Nqp;
        end

        %% Function to compute the evaluation function
        %
        %       Calling this function will call the subclass 'phi'
        %       method evaluated at some position vector, perform 
        %       an inner product with the weights, and return the 
        %       evaluated function as the combination of the
        %       Basis Functions in the set.
        %
        %   @param r Position vector (Rxd) where d is position
        %   dimensionality and R is number of positions
        %
        %   @return u Function evaluations at all positions (1xR)
        function u = u(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'BasisSet')}
                r (:,:) {mustBeNumeric}
            end
            % Sum of c*phi, where c (Nx1) and phi (NxR)
            u = obj.c'*obj.phi(r); 
        end

        %% Function to compute the evaluation of the 1st derivative
        %
        %       Calling this function will call the subclass 'dphi'
        %       method evaluated at some position vector, perform 
        %       an inner product with the weights, and return the 
        %       evaluated functions 1st derivative as the combination of
        %       the derivatives of the Basis Functions in the set.
        %
        %   @param r Position vector (Rxd) where d is position
        %   dimensionality and R is number of positions
        %
        %   @return du Function evaluations at all positions (1xR)
        function du = du(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'BasisSet')}
                r (:,:) {mustBeNumeric}
            end
            % Sum of c*dphi, where c (Nx1) and phi (NxR)
            du = obj.c'*obj.dphi(r); 
        end

        %% Function to compute the evaluation of the 2nd derivative
        %
        %       Calling this function will call the subclass 'd2phi'
        %       method evaluated at some position vector, perform 
        %       an inner product with the weights, and return the 
        %       evaluated functions 2nd derivative as the combination of
        %       the derivatives of the Basis Functions in the set.
        %
        %   @param r Position vector (Rxd) where d is position
        %   dimensionality and R is number of positions
        %
        %   @return d2u Function evaluations at all positions (1xR)
        function d2u = d2u(obj,r)
            arguments
                obj (1,1) {mustBeA(obj,'BasisSet')}
                r (:,:) {mustBeNumeric}
            end
            % Sum of c*d2phi, where c (Nx1) and phi (NxR)
            d2u = obj.c'*obj.d2phi(r); 
        end

    end

    methods (Abstract)
        phi(obj,r) % NxR matrix of N basis func at R pos
        dphi(obj,r) % NxR matrix of N derivatives of func at R pos
        d2phi(obj,r) % NxR matrix of N 2nd derivatives of func at R pos

        % NjxNi Overlap integral with phi (where Nj=size of rhs_obj and Ni=size of lhs_obj)
        integral_x_phi(rhs_obj, lhs_obj, lb, ub); 

        % NjxNi Overlap integral with dphi (where Nj=size of rhs_obj and Ni=size of lhs_obj)
        integral_x_dphi(rhs_obj, lhs_obj, lb, ub); 

        % NjxNi Overlap integral with dphi (where Nj=size of rhs_obj and Ni=size of lhs_obj)
        integral_x_d2phi(rhs_obj, lhs_obj, lb, ub);
    end
end