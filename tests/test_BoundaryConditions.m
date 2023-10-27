% Unit and integration tests for the BoundaryCondition class

%% Test 1 - Scalar Dirichlet BC 

vset = Variables();
u1 = Monomial1D(2);
u2 = Monomial1D(5);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;...
                            4;5;6;7;8]);


c_1 = vset.get_variable_coefficients('u1');
c_2 = vset.get_variable_coefficients('u2');

l1 = 1;
l2 = 2;

bc1 = DirichletBC(1);
bc2 = DirichletBC(0);

r11 = bc1.eval_dLdc(u1,c_1,1,l1);
r12 = bc1.eval_dLdl(u1,c_1,1,l1);

assert( norm(r11-[1;1]) < 1e-4 )
assert( abs(r12 + 2) < 1e-4 )

r21 = bc1.eval_dLdc(u2,c_2,0,l2);
r22 = bc1.eval_dLdl(u2,c_2,0,l2);

assert( norm(r21-[2;0;0;0;0]) < 1e-4 )
assert( abs(r22 - 3) < 1e-4 )
