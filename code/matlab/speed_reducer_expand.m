function x = speed_reducer_expand(xsFree, p, fixed)
% SPEED_REDUCER_EXPAND  Scaled free variables -> full x in display units.
%   fixed is [] or a struct with fields idx (indices) and val (display values).
%
% Author: AWD Labs
% Student ID: 52104479

xsFree = xsFree(:);
if isempty(fixed)
    x = p.lb + xsFree.*(p.ub - p.lb);
    return
end
free = setdiff(1:7, fixed.idx);
xs = zeros(7, 1, 'like', xsFree);
xs(free) = xsFree;
x = p.lb + xs.*(p.ub - p.lb);
x(fixed.idx) = fixed.val(:);
end
