function J = complex_step_jac(fun, x, h)
% COMPLEX_STEP_JAC  Jacobian of a real analytic function by complex step.
%   J(i,j) = Im(fun_i(x + i h e_j))/h, accurate to machine precision.
%
% Author: AWD Labs
% Student ID: 52104479

if nargin < 3, h = 1e-30; end
x = x(:);
n = numel(x);
f0 = fun(complex(x));
J = zeros(numel(f0), n);
for k = 1:n
    xc = complex(x);
    xc(k) = xc(k) + 1i*h;
    J(:, k) = imag(fun(xc))/h;
end
end
