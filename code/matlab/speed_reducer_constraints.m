function [c, ceq, gc, gceq] = speed_reducer_constraints(xs, p, fixed)
% SPEED_REDUCER_CONSTRAINTS  Normalised inequality constraints c <= 0 (11 rows).
%
% Author: AWD Labs
% Student ID: 52104479
%
% Stress, deflection and overall-size rows are value/limit - 1. The four
% geometric rows are linear in x, so they use a constant reference scale
% (smallest allowed length or width) to stay linear and convex.
% Row order: sigma_b, sigma_c, sigma_s1, sigma_s2, y1, y2, l1min, l2min,
%            b/m lower, b/m upper, overall size.

if nargin < 3, fixed = []; end
x = speed_reducer_expand(xs, p, fixed);
a = speed_reducer_core(x, p);
refL = p.lb(4)*p.cm;
refB = p.lb(1)*p.cm;
c = [a.sigma_b/p.sigma_b_max - 1;
     a.sigma_c/p.sigma_c_max - 1;
     a.sigma_s(1)/p.sigma_s_max - 1;
     a.sigma_s(2)/p.sigma_s_max - 1;
     a.y(1)/p.y_max - 1;
     a.y(2)/p.y_max - 1;
     (p.l1_min_factor*a.d(1) + p.l1_min_offset - a.l(1))/refL;
     (p.l2_min_factor*a.d(2) + p.l2_min_offset - a.l(2))/refL;
     (p.bm_min*a.m - a.b)/refB;
     (a.b - p.bm_max*a.m)/refB;
     a.overall/p.overall_max - 1];
ceq = [];
gceq = [];
if nargout > 2
    J = complex_step_jac(@(v) local_c(v, p, fixed), xs);
    gc = J.';          % fmincon wants n-by-m
end
end

function c = local_c(v, p, fixed)
c = speed_reducer_constraints(v, p, fixed);
end
