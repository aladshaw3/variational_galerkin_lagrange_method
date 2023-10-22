% Basic tests for overlap

%% Test 1 - Monomial1D overlap

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(2);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;-3;4;5]);

overlap = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u2.phi(r), -1, 1);

overlap_exact = integral_x_phi(vset.get_variable('u1'), u2, -1, 1);

overlap2 = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u1.phi(r), 0, 10);

overlap2_exact = vset.get_variable('u1').integral_x_phi(u1, 0, 10);

rel_error = norm(overlap - overlap_exact)/norm(overlap_exact);
rel_error2 = norm(overlap2 - overlap2_exact)/norm(overlap2_exact);

assert( rel_error < 1e-3 )
assert( rel_error2 < 1e-3 )

% Size of overlap:  Rows = size of j, Cols = size of i
%
%   overlap_ji = phi_j * phi_i

% Overlap returns an N1xN2 matrix

%% Test 2 - Gradient Overlaps (phi x dphi)

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(2);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;-3;4;5]);

overlap = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u2.dphi(r), -1, 1);

overlap_exact = integral_x_dphi(vset.get_variable('u1'), u2, -1, 1);

overlap2 = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u1.dphi(r), 0, 10);

overlap2_exact = vset.get_variable('u1').integral_x_dphi(u1, 0, 10);

rel_error = norm(overlap - overlap_exact)/norm(overlap_exact);
rel_error2 = norm(overlap2 - overlap2_exact)/norm(overlap2_exact);

assert( rel_error < 1e-3 )
assert( rel_error2 < 1e-3 )


%% Test 3 - Laplacian Overlaps (phi x d2phi)

vset = Variables();
u1 = Monomial1D(3);
u2 = Monomial1D(2);

vset.add_variable(u1);
vset.add_variable(u2);

vset.set_all_coefficients([1;-2;-3;4;5]);

overlap = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u2.d2phi(r), -1, 1);

overlap_exact = integral_x_d2phi(vset.get_variable('u1'), u2, -1, 1);

overlap2 = OverlapIntegrals_Numerical.phi_x_fun(vset.get_variable('u1'), @(r)u1.d2phi(r), 0, 10);

overlap2_exact = vset.get_variable('u1').integral_x_d2phi(u1, 0, 10);

error = norm(overlap - overlap_exact);
rel_error2 = norm(overlap2 - overlap2_exact)/norm(overlap2_exact);

assert( error < 1e-3 )
assert( rel_error2 < 1e-3 )