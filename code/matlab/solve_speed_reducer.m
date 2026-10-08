function r = solve_speed_reducer(method, xs0, p, fixed, record)
% SOLVE_SPEED_REDUCER  Run fmincon ('sqp' or 'interior-point') on the scaled problem.
%
% Author: AWD Labs
% Student ID: 52104479
%
%   r = solve_speed_reducer(method, xs0, p, fixed, record)
%   xs0    scaled start for the FREE variables (column)
%   fixed  [] or struct(idx, val) for fixed variables (display units)
%   record true to store per-iteration history in r.history = [f_cm3, maxViol]
%
% Options (tight tolerances, reported in the results):
%   OptimalityTolerance 1e-10, ConstraintTolerance 1e-10, StepTolerance 1e-12,
%   MaxIterations 1000, MaxFunctionEvaluations 1e5, user supplied gradients.

if nargin < 4, fixed = []; end
if nargin < 5, record = false; end
if ~any(strcmp(method, {'sqp','interior-point'}))
    error('solve_speed_reducer:badMethod', 'method must be "sqp" or "interior-point".');
end

nFree = numel(xs0);
history = [];
    function stop = outfun(~, optimValues, state)
        stop = false;
        if record && strcmp(state, 'iter')
            history(end+1, :) = [optimValues.fval*p.f_scale/p.cm^3, ...
                                 max(0, optimValues.constrviolation)]; %#ok<AGROW>
        end
    end

opts = optimoptions('fmincon', 'Algorithm', method, 'Display', 'off', ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
    'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10, ...
    'StepTolerance', 1e-12, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 1e5, ...
    'OutputFcn', @outfun);

objFun = @(v) speed_reducer_objective(v, p, fixed);
conFun = @(v) speed_reducer_constraints(v, p, fixed);

tic;
[xs, fval, exitflag, output, lambda] = fmincon(objFun, xs0(:), [], [], [], [], ...
    zeros(nFree,1), ones(nFree,1), conFun, opts);
r.time = toc;

c = conFun(xs);
r.method = method;
r.xs = xs;
r.x = speed_reducer_expand(xs, p, fixed);
r.f = fval*p.f_scale/p.cm^3;
r.viol = max(0, max(c));
r.feasible = r.viol <= 1e-6;
r.exitflag = exitflag;
r.ok = (exitflag > 0) && r.feasible;
r.nit = output.iterations;
r.nfev = output.funcCount;
r.firstorderopt = output.firstorderopt;
r.lambda = lambda;
r.message = output.message;
r.history = history;
end
