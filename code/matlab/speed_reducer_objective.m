function [f, g] = speed_reducer_objective(xs, p, fixed)
% SPEED_REDUCER_OBJECTIVE  Scaled volume (V / f_scale) and complex step gradient.
%
% Author: AWD Labs
% Student ID: 52104479

if nargin < 3, fixed = []; end
a = speed_reducer_core(speed_reducer_expand(xs, p, fixed), p);
f = a.V_total/p.f_scale;
if nargout > 1
    g = complex_step_jac(@(v) speed_reducer_objective(v, p, fixed), xs);
    g = g(:);
end
end
