%% Test 1 - Evaluation of Diffusion Residual

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(5);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;-3; ...
                            4;5;6;7;8]);


c_1 = vset.get_variable_coefficients('u1');
c_2 = vset.get_variable_coefficients('u2');

lb = 0;
ub = 1;
D = 0.5;
diff = Diffusion(lb,ub,D);

Res = diff.eval(u1,u2,c_2);
Res2 = diff.eval(u1,u1,c_1);
Res3 = diff.eval(u2,u1,c_1);

ol = -D*(vset.get_variable('u1').integral_x_d2phi(u2, 0, 1) * c_2);

er = norm(Res - ol)/norm(ol);

assert( er < 1e-3 )