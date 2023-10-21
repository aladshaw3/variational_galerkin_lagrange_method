% Unit tests for the Variables object

%% Test 1 - Basic Usage

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(2);

vset.add_variable(u1);
vset.add_variable(u2);

% Demonstrate that changes to Base u1 
% also changes u1 held in vset
u1.set_coeff([1;2;3]);
assert( norm(u1.get_coeff() - vset.get_variable('u1').get_coeff()) < 1e-6 )

c = vset.get_all_coefficients();
vset.set_all_coefficients([1;-2;-3;4;5]);

c = vset.get_all_coefficients();
c1 = u1.get_coeff();
assert( norm(c(1:3) - c1) < 1e-6 )

vset.set_variable_coefficients('u2',[-5;-4]);

c2 = vset.get_variable_coefficients('u2');
c2new = u2.get_coeff();
c = vset.get_all_coefficients();
assert( norm(c2new - c2) < 1e-6 )
