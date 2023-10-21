% Basic tests for overlap

%% Test 1 - Monomial1D overlap

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(2);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;-3;4;5]);

overlap = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u2.phi(r), -1, 1);

overlap_exact = OverlapIntegrals_Exact.phi_x_phi_Monomial1D_x_Monomial1D(vset.get_variable('u1'), u2, -1, 1);

overlap2 = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u1.phi(r), 0, 10);

overlap2_exact = OverlapIntegrals_Exact.phi_x_phi_Monomial1D_x_Monomial1D(vset.get_variable('u1'), u1, 0, 10);

rel_error = norm(overlap - overlap_exact)/norm(overlap_exact);
rel_error2 = norm(overlap2 - overlap2_exact)/norm(overlap2_exact);

assert( rel_error < 1e-3 )
assert( rel_error2 < 1e-3 )

% Size of overlap:  Rows = size of j, Cols = size of i
%
%   overlap_ji = phi_j * phi_i

% Overlap returns an N1xN2 matrix