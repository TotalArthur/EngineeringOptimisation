% POSTPROCESS_DESIGN  Constraint and bound tables, multipliers, Hessian, convexity.
%
% Author: AWD Labs
% Student ID: 52104479
%
% Run after run_optimisation. Reads results/matlab/multistart.mat and writes
% results/matlab/postprocess.mat plus csv tables. Multipliers and first order
% optimality come straight from the fmincon output of the best run.

clear; clc;
outDir = fullfile('..','..','results','matlab');
S = load(fullfile(outDir, 'multistart.mat'));
p = S.p;
methods = {'sqp', 'interior_point'};
activeTol = 1e-4;     % relative physical margin below this is "active"
boundTol = 1e-4;      % distance from a bound as a fraction of its range
nSeg = 5000;

post = struct();
for im = 1:2
    m = methods{im};
    b = S.best.(m);
    x = b.x;  xs = b.xs;
    a = speed_reducer_analysis(x, p);
    phys = physicalConstraints(x, a, p);
    margin = 100*(phys.limit - phys.value)./abs(phys.limit);
    active = (phys.limit - phys.value)./abs(phys.limit) < activeTol;
    lamC = b.lambda.ineqnonlin(:);

    % bounds table
    status = repmat({'free'}, 7, 1);
    status(xs < boundTol) = {'lower'};
    status(xs > 1 - boundTol) = {'upper'};

    % Hessian of the Lagrangian (scaled space) by central differences of
    % complex step gradients of L = f + lambda' c
    H = zeros(7);  h = 1e-6;
    for j = 1:7
        e = zeros(7,1);  e(j) = h;
        H(:,j) = (lagGrad(xs+e, p, lamC) - lagGrad(xs-e, p, lamC))/(2*h);
    end
    H = (H + H.')/2;
    eigH = sort(eig(H));

    % reduced Hessian on the null space of the active gradients
    Jc = complex_step_jac(@(v) speed_reducer_constraints(v, p, []), xs);
    A = Jc(active, :).';
    for i = find(xs < boundTol).',  e = zeros(7,1); e(i) = -1; A = [A e]; end %#ok<AGROW>
    for i = find(xs > 1-boundTol).', e = zeros(7,1); e(i) = 1;  A = [A e]; end %#ok<AGROW>
    if isempty(A), Z = eye(7); else, Z = null(A.'); end
    if isempty(Z), eigZ = []; else, eigZ = sort(eig(Z.'*H*Z)); end

    post.(m).x = x;  post.(m).V_cm3 = a.display.V_cm3;
    post.(m).constraints = phys;  post.(m).margin_pct = margin;  post.(m).active = active;
    post.(m).lambda_ineq = lamC;
    post.(m).lambda_lower = b.lambda.lower;  post.(m).lambda_upper = b.lambda.upper;
    post.(m).bound_status = status;
    post.(m).firstorderopt = b.firstorderopt;
    post.(m).eigH = eigH;  post.(m).eigHred = eigZ;
    post.(m).all_satisfied = all(phys.value <= phys.limit*(1 + 1e-9) + 1e-12);

    fprintf('\n==== %s: V = %.4f cm^3, first-order optimality %.2e ====\n', m, a.display.V_cm3, b.firstorderopt);
    for k = 1:numel(phys.name)
        fprintf('  %-34s %10.5g %-3s limit %10.5g margin %8.3f %%  %s  lam = %.4g\n', phys.name{k}, ...
            phys.value(k), phys.unit{k}, phys.limit(k), margin(k), ternary(active(k),'ACTIVE','inactive'), lamC(k));
    end
    for i = 1:7
        fprintf('  %-3s %9.5f %-3s [%g, %g] %s\n', p.var_names{i}, x(i), p.var_units{i}, p.lb(i), p.ub(i), status{i});
    end
    fprintf('  eig(H_L): %s\n', mat2str(eigH.', 4));
    fprintf('  eig(reduced H_L), dim %d: %s\n', size(Z,2), mat2str(eigZ.', 4));
    fprintf('  all original constraints satisfied: %d\n', post.(m).all_satisfied);

    T = table(phys.name, phys.value, phys.limit, phys.unit, margin, active, lamC, ...
        'VariableNames', {'constraint','value','limit','unit','margin_pct','active','multiplier'});
    writetable(T, fullfile(outDir, ['constraints_' m '.csv']));
    Tb = table(p.var_names.', x, p.lb, p.ub, status, 'VariableNames', {'variable','value','lb','ub','status'});
    writetable(Tb, fullfile(outDir, ['bounds_' m '.csv']));
end

%% sampled midpoint convexity test (objective and the 11 constraints)
rng(2);
viol = zeros(12,1);  worst = zeros(12,1);
for s = 1:nSeg
    xa = rand(7,1);  xb = rand(7,1);  xm = (xa + xb)/2;
    fa = [speed_reducer_objective(xa, p); speed_reducer_constraints(xa, p)];
    fb = [speed_reducer_objective(xb, p); speed_reducer_constraints(xb, p)];
    fm = [speed_reducer_objective(xm, p); speed_reducer_constraints(xm, p)];
    gap = fm - (fa + fb)/2;
    sc = max(1e-12, (abs(fa) + abs(fb))/2);
    viol = viol + (gap > 1e-9*sc);
    worst = max(worst, gap./sc);
end
names = [{'objective'}; {'gear bending stress'; 'gear contact stress'; 'shaft 1 combined stress'; ...
    'shaft 2 combined stress'; 'shaft 1 deflection'; 'shaft 2 deflection'; 'shaft 1 min length'; ...
    'shaft 2 min length'; 'b/m lower limit'; 'b/m upper limit'; 'overall size'}];
fprintf('\nMidpoint convexity test, %d segments\n', nSeg);
for k = 1:12
    fprintf('  %-26s violations %5d (%.1f %%) worst rel gap %.3g\n', names{k}, viol(k), 100*viol(k)/nSeg, worst(k));
end
post.convexity = table(names, viol, viol/nSeg, worst, 'VariableNames', {'function','n_violations','fraction','worst_rel_gap'});
writetable(post.convexity, fullfile(outDir, 'convexity_samples.csv'));
save(fullfile(outDir, 'postprocess.mat'), 'post');

%% local functions
function g = lagGrad(xs, p, lam)
J = complex_step_jac(@(v) speed_reducer_constraints(v, p, []), xs);
g = speed_reducer_objective_grad(xs, p) + J.'*lam;
end

function g = speed_reducer_objective_grad(xs, p)
[~, g] = speed_reducer_objective(xs, p, []);
end

function phys = physicalConstraints(x, a, p)
% Original, unnormalised constraints written out directly (independent of c).
d = a.display;  cm = p.cm;
phys.name = {'gear bending stress'; 'gear contact stress'; 'shaft 1 combined stress'; ...
    'shaft 2 combined stress'; 'shaft 1 deflection'; 'shaft 2 deflection'; ...
    'shaft 1 min length (1.5 d1 + 1.9)'; 'shaft 2 min length (1.1 d2 + 1.9)'; ...
    'b/m lower limit'; 'b/m upper limit'; 'overall size m z (1+u)'};
phys.value = [d.sigma_b_MPa; d.sigma_c_MPa; d.sigma_s_MPa(1); d.sigma_s_MPa(2); d.y_mm(1); d.y_mm(2); ...
    p.l1_min_factor*x(6) + p.l1_min_offset/cm; p.l2_min_factor*x(7) + p.l2_min_offset/cm; ...
    p.bm_min; a.b_over_m; a.overall/cm];
phys.limit = [p.sigma_b_max/p.MPa; p.sigma_c_max/p.MPa; p.sigma_s_max/p.MPa; p.sigma_s_max/p.MPa; ...
    p.y_max/p.mm; p.y_max/p.mm; x(4); x(5); a.b_over_m; p.bm_max; p.overall_max/cm];
phys.unit = {'MPa'; 'MPa'; 'MPa'; 'MPa'; 'mm'; 'mm'; 'cm'; 'cm'; '-'; '-'; 'cm'};
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
