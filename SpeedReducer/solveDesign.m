function out = solveDesign(algorithm, x0, params, zFixed)
% solveDesign  Run fmincon once and collect the results in a struct.
% Author: AWD Labs
%
% INPUTS
%   algorithm  'sqp' or 'interior-point'
%   x0         starting point, brief units [b m z l1 l2 d1 d2]
%   params     parameter struct from main.m
%   zFixed     (optional) if given, z is held at this value and only the
%              other six variables are optimised (used for the rounded-z run)
%
% OUTPUT (struct)
%   x, V (cm^3), exitflag, iterations, funcCount, time (s),
%   maxViol (largest constraint value c, feasible if <= params.feasTol)
%
% Both algorithms get exactly the same options, objective and constraints,
% so the comparison is fair.

if nargin < 4
    zFixed = NaN;
end

opts = optimoptions('fmincon', ...
    'Algorithm',               algorithm, ...
    'Display',                 'off', ...
    'OptimalityTolerance',     params.optTol, ...
    'ConstraintTolerance',     params.conTol, ...
    'StepTolerance',           params.stepTol, ...
    'MaxIterations',           params.maxIter, ...
    'MaxFunctionEvaluations',  params.maxEval);

if isnan(zFixed)
    % all seven variables free
    expand = @(v) v;
    start  = x0;
    lb     = params.lb;
    ub     = params.ub;
else
    % z fixed: the solver only sees [b m l1 l2 d1 d2] and we put z back in
    keep   = [1 2 4 5 6 7];
    expand = @(v) [v(1:2), zFixed, v(3:6)];
    start  = x0(keep);
    lb     = params.lb(keep);
    ub     = params.ub(keep);
end

fun    = @(v) objectiveFun(expand(v), params);
nonlin = @(v) constraintFun(expand(v), params);

tStart = tic;
[v, fval, exitflag, output] = fmincon(fun, start, [], [], [], [], lb, ub, nonlin, opts);
runTime = toc(tStart);

out.x          = expand(v);
out.V          = fval * params.objScale;     % back to cm^3
out.exitflag   = exitflag;
out.iterations = output.iterations;
out.funcCount  = output.funcCount;
out.time       = runTime;
out.maxViol    = max(constraintFun(out.x, params));
end
