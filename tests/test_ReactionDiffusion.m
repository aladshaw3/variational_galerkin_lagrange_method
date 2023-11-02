% Now have the minimum of what is technically required to build and solve
% a Reaction diffusion equation...

% Solve:  -D * d2u/dx2 + k*u = 0   with u(0) = uo and du(L) = 0
%
%       Basis set: Monomial1D, with basis size 3 and 2 BCs

vset = Variables();
u = Monomial1D(3);
vset.add_variable(u);

vset.set_all_coefficients([0;0;0]);
c = vset.get_variable_coefficients('u');
% temporary storage of Lagrangian multipliers
l = [0;0];

D = 0.5;
k = 2;
L = 1;
uo = 1;
diff = Diffusion(0,L,D);
react = FirstOrderReaction(0,L,k);

bc1 = NeumannBC(0);
bc1_r = L;
bc2 = DirichletBC(uo);
bc2_r = 0;

x_span = 0:(L/10):L;
u_exact = uo*(exp((2*L-x_span)*sqrt(k/D))+exp(x_span*sqrt(k/D)))./(1+exp(2*L*sqrt(k/D)));

funop = @(x) diffreact(x,vset,diff,react,bc1,bc1_r,bc2,bc2_r);

% Alternative using fminsearch 
x = fminsearch(funop,[0;0;0;0;0]);

u_approx = u.u(x_span);

% Assert less than 5% error
assert( norm(u_exact-u_approx)/norm(u_exact) < 0.05 )

% Helper function for residuals
function Ax = diffreact(x,vset,diff,react,bc1,bc1_r,bc2,bc2_r)
    c1 = x(1);
    c2 = x(2);
    c3 = x(3);
    vset.set_all_coefficients([c1;c2;c3]);
    l1 = x(4);
    l2 = x(5);

    u = vset.get_variable('u');

    Ax = zeros(5,1);
    Ax(1:3,1) = diff.eval(u,u,[c1;c2;c3]) + ...
                    react.eval(u,u,[c1;c2;c3]) + ...
                    bc1.eval_dLdc(u,[c1;c2;c3],bc1_r,l1) + ...
                    bc2.eval_dLdc(u,[c1;c2;c3],bc2_r,l2);
    Ax(4,1) = bc1.eval_dLdl(u,[c1;c2;c3],bc1_r,l1);
    Ax(5,1) = bc2.eval_dLdl(u,[c1;c2;c3],bc2_r,l2);

    % Force a norm calculation (for fminsearch)
    Ax=norm(Ax);
end