% Now have the minimum of what is technically required to build and solve
% an steady-state Laplaces equation in 1D. So, let's test it out...

% Solve:  d2u/dx2 = 0   with u(0) = 1 and u(1) = 0
%
%       Basis set: Monomial1D, with basis size 2 and 2 BCs

vset = Variables();
u = Monomial1D(2);
vset.add_variable(u);

vset.set_all_coefficients([0;0]);
c = vset.get_variable_coefficients('u');
% temporary storage of Lagrangian multipliers
l = [0;0];

diff = Diffusion(0,1,1);
bc1 = DirichletBC(1);
bc1_r = 0;
bc2 = DirichletBC(0);
bc2_r = 1;

Ax = lindiff([c;l], vset, diff, bc1,bc1_r,bc2,bc2_r);

% Only valid for vector form with gmres
%assert( norm(Ax - [0;0;-1;0]) < 1e-4 )

linop = @(x) lindiff(x,vset,diff,bc1,bc1_r,bc2,bc2_r);
%b = [0;0;0;1e-16];  % NOTE: Cannot give b = zeros for gmres, so this 
%                    % is used to force an evaluation
%x = gmres(linop,b);

% Alternative using fminsearch 
x = fminsearch(linop,[0;0;0;0]);

assert( norm(x - [1;-1;0;0]) < 1e-4 )

dist = 0:0.1:1;
usol = u.u(dist);

% Assert solution is correct 
assert( norm(usol - [1,.9,.8,.7,.6,.5,.4,.3,.2,.1,0]) < 1e-4 )



% Helper function for linear operator 
function Ax = lindiff(x,vset,diff,bc1,bc1_r,bc2,bc2_r)
    c1 = x(1);
    c2 = x(2);
    vset.set_all_coefficients([c1;c2]);
    l1 = x(3);
    l2 = x(4);

    u = vset.get_variable('u');

    Ax = zeros(4,1);
    Ax(1:2,1) = diff.eval(u,u,[c1;c2]) + ...
                    bc1.eval_dLdc(u,[c1;c2],bc1_r,l1) + ...
                    bc2.eval_dLdc(u,[c1;c2],bc2_r,l2);
    Ax(3,1) = bc1.eval_dLdl(u,[c1;c2],bc1_r,l1);
    Ax(4,1) = bc2.eval_dLdl(u,[c1;c2],bc2_r,l2);

    % Force a norm calculation (for fminsearch)
    Ax=norm(Ax);
end
