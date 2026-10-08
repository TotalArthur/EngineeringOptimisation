function [c, ceq] = constraintFun(x, params)
% constraintFun  Nonlinear constraints for fmincon, in the form c <= 0.
% Author: AWD Labs
%
% Each constraint is written as a utilisation = demand / capacity, so
%   utilisation <= 1  means the design is OK,
%   utilisation  = 1  means the constraint is exactly on its limit.
% fmincon gets c = utilisation - 1, so c <= 0 is feasible. Doing it this way
% puts every constraint on the same scale (stresses are around 1e8 Pa, but
% ratios are around 1), which helps the solver a lot.
%
% For a maximum limit (stress <= limit):  utilisation = value / limit
% For a minimum limit (value >= limit):   utilisation = limit / value

r = speedReducerAnalysis(x, params);

utilisation = r.conValue ./ r.conLimit;
isMin = logical(params.conIsMin);
utilisation(isMin) = r.conLimit(isMin) ./ r.conValue(isMin);

c   = utilisation - 1;
ceq = [];      % no equality constraints in this problem
end
