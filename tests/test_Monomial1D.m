% Unit tests for Monomial1D basis sets

%% Test 1 - Basic usage

% Create class and set parameters
obj = Monomial1D(3);

Nqp = obj.get_Nqp();
assert( norm([1;3;5] - Nqp) < 1e-6 )

c = [0; 0.5; 2];
obj.set_coeff(c);
assert( norm(c - obj.get_coeff()) < 1e-6 )

% Evaluation of basis functions at specific points
x = (1:5)';
phi = obj.phi(x);

% Calculate manually for confirmation
phi_test = zeros(length(c),length(x));
for i=1:length(c)
    for j=1:5
        phi_test(i,j) = j^(i-1);
    end
end
assert( norm(phi_test - phi) < 1e-6 )

% Evaluate function u = sum(c*phi) at each position r
u = obj.u(x);
% Result should be u = 0.5*x + 2*x^2
u_test = zeros(1,length(x));
for i=1:5
    u_test(i) = 0.5*x(i)+2*x(i)^2;
end
assert( norm(u_test - u) < 1e-6 )


% Estimate 1st derivatives
dphi = obj.dphi(x);
dphi_test = zeros(length(c),length(x));
for i=1:length(c)
    for j=1:5
        dphi_test(i,j) = (i-1)*j^(i-2);
    end
end
assert( norm(dphi_test - dphi) < 1e-6 )


% Evaluate function du = sum(c*dphi) at each position r
du = obj.du(x);
% Result should be du = 0.5 + 4*x
du_test = zeros(1,length(x));
for i=1:5
    du_test(i) = 0.5+4*x(i);
end
assert( norm(du_test - du) < 1e-6 )


% Estimate 2nd derivatives
d2phi = obj.d2phi(x);
d2phi_test = zeros(length(c),length(x));
for i=1:length(c)
    for j=1:5
        d2phi_test(i,j) = (i-2)*(i-1)*j^(i-3);
    end
end
assert( norm(d2phi_test - d2phi) < 1e-6 )


% Evaluate function du = sum(c*d2phi) at each position r
d2u = obj.d2u(x);
% Result should be d2u = 4
d2u_test = zeros(1,length(x));
for i=1:5
    d2u_test(i) = 4;
end
assert( norm(d2u_test - d2u) < 1e-6 )