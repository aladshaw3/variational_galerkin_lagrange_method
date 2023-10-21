classdef Variables < handle

    properties (Access = private)
        var_map         % Matlab struct of variable names and basis sets
        coeffs          % Nx1 vector of stacked coefficients
        cs_idx          % vector of starting indices for coefficients
    end

    methods (Access = public)
        %% Constructor 
        function obj = Variables()
            obj.var_map = struct();
            obj.coeffs = [];
            obj.cs_idx = [1];
        end

        %% Function to add a basis set to the variable set
        function [] = add_variable(obj, basis_set)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
                basis_set (1,1) {mustBeA(basis_set,'BasisSet')}
            end
            obj.var_map.(inputname(2)) = basis_set;
            obj.coeffs = [obj.coeffs; basis_set.get_coeff()];
            ids = length(obj.coeffs) + 1;
            obj.cs_idx = [obj.cs_idx; ids];
        end

        %% Function to return the basis set in struct by name
        function basis = get_variable(obj, name)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
                name (1,:) {mustBeTextScalar}
            end
            basis = obj.var_map.(name);
        end

        %% Function to return current set of coefficients
        function c = get_all_coefficients(obj)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
            end
            c = obj.coeffs();
            fnames = fieldnames(obj.var_map);
            idx=1;
            for i=1:numel(fnames)
                N = length(obj.var_map.(fnames{i}).get_coeff());
                c(idx:idx+N-1,1) = obj.var_map.(fnames{i}).get_coeff();
                idx=idx+N;
            end
        end

        %% Function to get the coefficients of a specific variable
        function c = get_variable_coefficients(obj, name)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
                name (1,:) {mustBeTextScalar}
            end
            c = obj.var_map.(name).get_coeff();
        end

        %% Function to set the coefficients for all variables
        function [] = set_all_coefficients(obj, c)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
                c (:,1) {mustBeNumeric}
            end
            if (length(c) ~= length(obj.coeffs))
                ME = MException('Variables:invalidSize',...
                    '\nInput vector must be of same size as "obj.coeffs", N=%g',string(length(obj.coeffs)));
                throw(ME)
            end
            obj.coeffs = c;
            fnames = fieldnames(obj.var_map);
            idx=1;
            for i=1:numel(fnames)
                N = length(obj.var_map.(fnames{i}).get_coeff());
                obj.var_map.(fnames{i}).set_coeff(c(idx:idx+N-1,1));
                idx=idx+N;
            end
        end

        %% Function to set coefficients of specific variable
        function [] = set_variable_coefficients(obj, name, c)
            arguments
                obj (1,1) {mustBeA(obj,'Variables')}
                name (1,:) {mustBeTextScalar}
                c (:,1) {mustBeNumeric}
            end
            if (length(c) ~= length(obj.var_map.(name).get_coeff()))
                ME = MException('Variables:invalidSize',...
                    '\nInput vector must be of same size as "obj.coeffs", N=%g',string(length(obj.coeffs)));
                throw(ME)
            end
            obj.var_map.(name).set_coeff(c);
            fnames = fieldnames(obj.var_map);
            idx = find(ismember(fnames, name));
            num_vals = obj.cs_idx(idx+1) - obj.cs_idx(idx);
            obj.coeffs(obj.cs_idx(idx):num_vals+obj.cs_idx(idx)-1) = c;
        end
    end
end