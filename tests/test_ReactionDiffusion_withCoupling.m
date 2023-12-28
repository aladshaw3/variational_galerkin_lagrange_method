% Now have the minimum of what is technically required to build and solve
% a Reaction diffusion equation...

% Solve:        -D * d2u1/dx2 + k*u1 = 0   with u1(0) = u1o and du1(L) = 0
% coupled with: -D * d2u2/dx2 - k*u1 = 0   with u2(0) = u2o and du2(L) = 0
%
%   Basis set: Monomial1D, with basis size 3 (x2 vars) and 2 BCs (x2 vars)

vset = Variables();
u1 = Monomial1D(4);
u2 = Monomial1D(4);
vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([0;0;0;0;0;0;0;0]);

% temporary storage of Lagrangian multipliers
l = [0;0;0;0];

D = 0.5;
k = 2;
L = 1;
u1o = 1;
u2o = 0;
diff1 = Diffusion(0,L,D);
diff2 = Diffusion(0,L,D);
react1 = FirstOrderReaction(0,L,k);
react2 = FirstOrderReaction(0,L,-k);

bc1_1 = NeumannBC(0);
bc1_1r = L;
bc2_1 = DirichletBC(u1o);
bc2_1r = 0;

bc1_2 = NeumannBC(0);
bc1_2r = L;
bc2_2 = DirichletBC(u2o);
bc2_2r = 0;

x_span = 0:(L/10):L;
u1_exact = u1o*(exp((2*L-x_span)*sqrt(k/D))+exp(x_span*sqrt(k/D)))./(1+exp(2*L*sqrt(k/D)));

funop = @(x) coupleddiffreact(x,vset,diff1,react1,bc1_1,bc1_1r,bc2_1,bc2_1r,...
                                    diff2,react2,bc1_2,bc1_2r,bc2_2,bc2_2r);

% Alternative using fminsearch 
options = optimset('MaxFunEvals',100000,'Display','iter','MaxIter',10000,'TolFun',1e-14,'TolX',1e-14);
x = fminsearch(funop,[0;0;0;0;0;0;0;0;0;0;0;0],options);

u1_approx = u1.u(x_span);
u2_approx = u2.u(x_span);
u_sum = u1_approx+u2_approx;

% Assert the expected result
assert( norm(u1_exact-u1_approx)/norm(u1_exact) < 0.05 )

% Assert conservation of mass 
assert( norm(u_sum - 1) < 0.05 )

% Helper function for residuals
function Ax = coupleddiffreact(x,vset,diff1,react1,bc1_1,bc1_1r,bc2_1,bc2_1r,...
                                    diff2,react2,bc1_2,bc1_2r,bc2_2,bc2_2r)
    full_c = x(1:8);

    vset.set_all_coefficients(full_c);
    c1_set = vset.get_variable_coefficients('u1');
    c2_set = vset.get_variable_coefficients('u2');


    l1 = x(9);
    l2 = x(10);
    l3 = x(11);
    l4 = x(12);

    u1 = vset.get_variable('u1');
    u2 = vset.get_variable('u2');

    Ax = zeros(12,1);
    Ax(1:4,1) = diff1.eval(u1,u1,c1_set) + ...
                    react1.eval(u1,u1,c1_set) + ...
                    bc1_1.eval_dLdc(u1,c1_set,bc1_1r,l1) + ...
                    bc2_1.eval_dLdc(u1,c1_set,bc2_1r,l2);

    Ax(5:8,1) = diff2.eval(u2,u2,c2_set) + ...
                    react2.eval(u2,u1,c1_set) + ...
                    bc1_2.eval_dLdc(u2,c2_set,bc1_2r,l3) + ...
                    bc2_2.eval_dLdc(u2,c2_set,bc2_2r,l4);

    Ax(9,1) = bc1_1.eval_dLdl(u1,c1_set,bc1_1r,l1);
    Ax(10,1) = bc2_1.eval_dLdl(u1,c1_set,bc2_1r,l2);

    Ax(11,1) = bc1_2.eval_dLdl(u2,c2_set,bc1_2r,l3);
    Ax(12,1) = bc2_2.eval_dLdl(u2,c2_set,bc2_2r,l4);

    % Force a norm calculation (for fminsearch)
    Ax=norm(Ax);
end