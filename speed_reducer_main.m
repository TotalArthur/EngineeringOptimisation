%% SPEED REDUCER OPTIMISATION (EG503X Assignment 1)
% Author: AWD Labs
% Student ID: 52104479
%
% Open this file in MATLAB and press Run. It does everything:
%   1. validates the analysis against the brief's validation case
%   2. solves the problem with fmincon 'sqp' and 'interior-point'
%   3. multistart from random points (rng seed 1) for both methods
%   4. integer teeth study (z = 17 to 28)
%   5. post processing of the best designs (constraints, bounds, multipliers,
%      Hessian of the Lagrangian, convexity test, sanity check)
%   6. reference case (u = 3, T = 1000 Nm) and comparison with Golinski (1970)
%   7. figures (PNG 300 dpi and PDF) and csv/mat results
%
% Needs the Optimization Toolbox (fmincon). Takes a few minutes to run.
% All the functions it uses are at the bottom of this file.
%
% UNITS
%   The design vector x is always in these units:
%       x = [b (cm); m (mm); z (-); l1 (cm); l2 (cm); d1 (cm); d2 (cm)]
%   The analysis function converts to SI (m, Pa, N m, m^3) straight away.
%   Results are converted back (cm^3, MPa, mm) only for display.

clear; clc; close all;

if exist('fmincon', 'file') == 0
    error('This script needs the Optimization Toolbox (fmincon was not found).');
end

%% 0. SETTINGS AND PARAMETERS  (change things here)

studentID = '52104479';
nStarts    = 100;    % random starts per method (continuous problem)
nStartsZ   = 20;     % random starts per method for each fixed z
nStartsRef = 50;     % random starts per method for the reference case
seed       = 1;      % random seed so results can be reproduced
methods    = {'sqp', 'interior-point'};
makeFigures = true;

% folders for outputs (next to this script)
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
resDir = fullfile(here, 'results');
figDir = fullfile(here, 'figures');
if ~exist(resDir, 'dir'), mkdir(resDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end
diary off; diary(fullfile(resDir, 'run_log.txt'));

% loading from the student ID: last digit n, second last digit m
n = str2double(studentID(end));
m = str2double(studentID(end-1));
p.T = 1500 + 100*n;          % input torque (Nm)
p.u = 2.5 + 0.1*m;           % speed ratio (-)

% material and gear constants
p.q  = 2.54;                 % tooth form factor
p.kv = 2.1;                  % contact stress concentration factor
p.E  = 200e9;                % Young's modulus (Pa)
p.nu = 0.3;                  % Poisson's ratio
% Elastic coefficient. The brief prints (1 + nu^2) but only (1 - nu^2)
% reproduces the validation value of 532 MPa (see section 1).
p.cp = p.E/(2*pi*(1 - p.nu^2));

% limits (SI)
p.sigb_max = 650e6;          % gear bending stress (Pa)
p.sigc_max = 800e6;          % gear contact stress (Pa)
p.sigs_max = 550e6;          % shaft combined stress (Pa)
p.y_max    = 0.075e-3;       % shaft deflection (m)

% bounds, in the units of x:  b    m    z   l1   l2   d1   d2
p.lb = [2.6; 7.0; 17; 7.3; 7.3; 2.8; 5.0];
p.ub = [4.4; 8.0; 28; 8.3; 8.3; 3.9; 5.5];
varNames = {'b', 'm', 'z', 'l1', 'l2', 'd1', 'd2'};
varUnits = {'cm', 'mm', '-', 'cm', 'cm', 'cm', 'cm'};

% gear cross section (Figure 2 of the brief) and geometric requirements
p.dst_off = 2.4;             % d_st = m (z - 2.4)
p.dw_off  = 4;               % d_w  = d_st - 4 m
p.dp_fac  = 2.4;             % d_p  = 2.4 d
p.web     = 1/3;             % web thickness = b/3
p.hub_len = 2;               % hub length = 2 d
p.l1_k = 1.5;  p.l1_c = 1.9; % 1.5 d1 + 1.9 <= l1   (1.9 in cm)
p.l2_k = 1.1;  p.l2_c = 1.9; % 1.1 d2 + 1.9 <= l2   (1.9 in cm)
p.bm_min = 5;  p.bm_max = 12;% 5 <= b/m <= 12
p.size_max = 160;            % m z (1+u) <= 160 cm

% unit conversions
p.cm = 1e-2;  p.mm = 1e-3;

% names for the 11 normalised constraints (same order as in constraints())
conNames = {'gear bending stress'; 'gear contact stress'; 'shaft 1 combined stress'; ...
    'shaft 2 combined stress'; 'shaft 1 deflection'; 'shaft 2 deflection'; ...
    'shaft 1 min length'; 'shaft 2 min length'; 'b/m lower limit'; ...
    'b/m upper limit'; 'overall size'};

fprintf('Student ID %s: T = %g Nm, u = %g\n\n', studentID, p.T, p.u);

%% 1. VALIDATION
% Validation data from the brief. Run at T = 1000 Nm (as printed) and at
% T = 2400 Nm to see which one reproduces the targets.

xVal = [3.5; 7.0; 22; 7.4; 7.8; 3.5; 5.2];     % b m z l1 l2 d1 d2
valNames = {'Volume (cm^3)', 'sigma_b (MPa)', 'sigma_c (MPa)', 'sigma_s1 (MPa)', ...
            'sigma_s2 (MPa)', 'y1 (mm)', 'y2 (mm)'};
valTarget = [4147; 323; 532; 512; 454; 0.0179; 0.0043];

valComp = zeros(7, 2);
torques = [1000, 2400];
for k = 1:2
    pv = p;  pv.T = torques(k);  pv.u = 3;     % validation case has u = 3
    r = analysis(xVal, pv);
    valComp(:, k) = [r.V/p.cm^3; r.sigb/1e6; r.sigc/1e6; r.sigs(1)/1e6; ...
                     r.sigs(2)/1e6; r.y(1)/p.mm; r.y(2)/p.mm];
end
valErr = 100*(valComp - valTarget)./valTarget;

fprintf('=== 1. VALIDATION ===\n');
fprintf('%-16s %10s | %10s %8s | %10s %8s\n', 'quantity', 'target', 'T=1000', 'err %', 'T=2400', 'err %');
for i = 1:7
    fprintf('%-16s %10.4g | %10.4g %8.2f | %10.4g %8.2f\n', valNames{i}, valTarget(i), ...
        valComp(i,1), valErr(i,1), valComp(i,2), valErr(i,2));
end

% which reading of cp matches the 532 MPa target (T = 2400 Nm)?
pv = p;  pv.T = 2400;  pv.u = 3;
pv.cp = pv.E/(2*pi*(1 - pv.nu^2));   rA = analysis(xVal, pv);
pv.cp = pv.E/(2*pi*(1 + pv.nu^2));   rB = analysis(xVal, pv);
fprintf('\nsigma_c with cp = E/(2 pi (1 - nu^2)): %.2f MPa  (target 532)\n', rA.sigc/1e6);
fprintf('sigma_c with cp = E/(2 pi (1 + nu^2)): %.2f MPa  (as printed in the brief)\n', rB.sigc/1e6);
fprintf('T = 2400 Nm matches every stress and deflection, T = 1000 Nm does not.\n');
fprintf('The volume is 4.5 %% below the target at both torques (unresolved, see summary).\n\n');

fid = fopen(fullfile(resDir, 'validation.csv'), 'w');
fprintf(fid, 'quantity,target,T1000,T1000_err_pct,T2400,T2400_err_pct\n');
for i = 1:7
    fprintf(fid, '%s,%g,%g,%g,%g,%g\n', valNames{i}, valTarget(i), valComp(i,1), valErr(i,1), valComp(i,2), valErr(i,2));
end
fclose(fid);

%% 2. FIRST SOLVE FROM THE CENTRE OF THE BOUNDS (both methods)
% fmincon with finite difference gradients. Tight tolerances are set inside
% run_solver (see the bottom of the file).

obj = @(x) objective(x, p);
con = @(x) constraints(x, p);
x0 = (p.lb + p.ub)/2;

fprintf('=== 2. SINGLE SOLVE FROM THE CENTRE OF THE BOUNDS ===\n');
for k = 1:2
    r = run_solver(methods{k}, obj, con, x0, p.lb, p.ub, 1000);
    fprintf('%-15s V = %.4f cm^3, exitflag %d, %d iterations, %d function evaluations, max violation %.1e\n', ...
        methods{k}, r.f, r.flag, r.iter, r.nfev, r.viol);
end
fprintf('\n');

%% 3. MULTISTART (continuous z)
% Random starting points inside the bounds. The problem is non-convex so
% one start can never prove a global optimum; many starts give evidence.

fprintf('=== 3. MULTISTART, %d starts per method (rng seed %d) ===\n', nStarts, seed);
rng(seed);
starts = zeros(7, nStarts);
for k = 1:nStarts
    starts(:, k) = p.lb + rand(7, 1).*(p.ub - p.lb);
end

runs = cell(1, 2);        % runs{method} is a struct array of all runs
best = cell(1, 2);        % best feasible run for each method
for k = 1:2
    for i = 1:nStarts
        runs{k}(i) = run_solver(methods{k}, obj, con, starts(:, i), p.lb, p.ub, 1000);
    end
    best{k} = best_run(runs{k});
    print_stats(methods{k}, runs{k});
    save_runs(fullfile(resDir, ['multistart_' strrep(methods{k}, '-', '_') '.csv']), runs{k});
end
fprintf('\nBest design, SQP:            %s\n', mat2str(best{1}.x.', 6));
fprintf('Best design, interior point: %s\n\n', mat2str(best{2}.x.', 6));

%% 4. INTEGER TEETH STUDY
% Fix z at each integer, re-optimise the other six variables (multistart).

fprintf('=== 4. INTEGER TEETH, z = 17 to 28 ===\n');
zList = 17:28;
Vz = nan(numel(zList), 2);
rng(seed);
lb6 = p.lb([1 2 4 5 6 7]);  ub6 = p.ub([1 2 4 5 6 7]);
for iz = 1:numel(zList)
    z = zList(iz);
    objZ = @(v) objective(insert_z(v, z), p);
    conZ = @(v) constraints(insert_z(v, z), p);
    for k = 1:2
        bestV = inf;
        for i = 1:nStartsZ
            v0 = lb6 + rand(6, 1).*(ub6 - lb6);
            r = run_solver(methods{k}, objZ, conZ, v0, lb6, ub6, 1000);
            if r.ok && r.f < bestV, bestV = r.f; end
        end
        if isfinite(bestV), Vz(iz, k) = bestV; end
    end
    fprintf('z = %2d:  SQP %9.3f cm^3   interior point %9.3f cm^3\n', z, Vz(iz,1), Vz(iz,2));
end
[Vint, iBest] = min(Vz(:, 1));
Vcont = best{1}.f;
fprintf('\nBest integer design: z = %d, V = %.3f cm^3\n', zList(iBest), Vint);
fprintf('Volume penalty of rounding vs continuous optimum: %.3f cm^3 (%.3f %%)\n\n', ...
    Vint - Vcont, 100*(Vint - Vcont)/Vcont);

fid = fopen(fullfile(resDir, 'integer_z.csv'), 'w');
fprintf(fid, 'z,V_sqp_cm3,V_interior_point_cm3\n');
for iz = 1:numel(zList)
    fprintf(fid, '%d,%.6f,%.6f\n', zList(iz), Vz(iz,1), Vz(iz,2));
end
fclose(fid);

%% 5. POST PROCESSING OF THE BEST DESIGNS
% Constraint table, bound table, multipliers, first order optimality,
% Hessian of the Lagrangian, and a check of the original constraints.

fprintf('=== 5. POST PROCESSING ===\n');
activeTol = 1e-4;    % active if the margin is below 0.01 % of the limit
boundTol  = 1e-4;    % at a bound if within 0.01 % of the variable range
margins = zeros(11, 2);
for k = 1:2
    b = best{k};
    x = b.x;
    [val, lim, unit] = physical_constraints(x, p);
    margin = 100*(lim - val)./abs(lim);          % % of limit still available
    isActive = (lim - val)./abs(lim) < activeTol;
    margins(:, k) = margin;

    fprintf('\n--- %s: V = %.4f cm^3 ---\n', methods{k}, b.f);
    fprintf('%-26s %10s %10s %5s %9s %-9s %10s\n', 'constraint', 'value', 'limit', 'unit', 'margin %', 'status', 'multiplier');
    for i = 1:11
        if isActive(i), st = 'ACTIVE'; else, st = 'inactive'; end
        fprintf('%-26s %10.5g %10.5g %5s %9.3f %-9s %10.4g\n', conNames{i}, val(i), lim(i), unit{i}, ...
            margin(i), st, b.lamC(i));
    end
    fprintf('\n%-4s %10s %5s %6s %6s %-8s\n', 'var', 'value', 'unit', 'lower', 'upper', 'status');
    atLower = (x - p.lb)./(p.ub - p.lb) < boundTol;
    atUpper = (p.ub - x)./(p.ub - p.lb) < boundTol;
    for i = 1:7
        if atLower(i), st = 'lower'; elseif atUpper(i), st = 'upper'; else, st = 'free'; end
        fprintf('%-4s %10.5f %5s %6.1f %6.1f %-8s\n', varNames{i}, x(i), varUnits{i}, p.lb(i), p.ub(i), st);
    end

    % sanity check: all original (unnormalised) constraints and bounds
    % (allowed to exceed a limit by 1e-6 of itself, the solver's feasibility tolerance)
    okCon = all(val <= lim*(1 + 1e-6));
    okBnd = all(x >= p.lb - 1e-6) && all(x <= p.ub + 1e-6);
    fprintf('\nAll original constraints satisfied: %d, all bounds satisfied: %d\n', okCon, okBnd);
    fprintf('Largest constraint overshoot: %.2e of the limit\n', max(0, max((val - lim)./abs(lim))));
    fprintf('First order optimality (fmincon output): %.3e\n', b.firstorder);
    % bound multipliers are per unit of each variable (cm, mm or teeth), so
    % compare them with care. Constraint multipliers above are dimensionless.
    fprintf('Lagrange multipliers on the lower bounds [b m z l1 l2 d1 d2]: %s\n', mat2str(b.lamLo.', 4));
    fprintf('Lagrange multipliers on the upper bounds [b m z l1 l2 d1 d2]: %s\n', mat2str(b.lamUp.', 4));

    % Hessian of the Lagrangian L = f + lambda'*c, worked out numerically.
    % Done in scaled variables xs in [0,1] so the eigenvalues are comparable.
    toX = @(xs) p.lb + xs.*(p.ub - p.lb);
    L   = @(xs) objective(toX(xs), p) + b.lamC.'*constraints(toX(xs), p);
    xs  = (x - p.lb)./(p.ub - p.lb);
    H   = num_hessian(L, xs, 1e-4);
    fprintf('Eigenvalues of the Hessian of the Lagrangian: %s\n', mat2str(sort(eig(H)).', 4));

    % reduced Hessian: only directions that keep the active constraints active
    cs = @(xs) constraints(toX(xs), p);
    J  = num_grad(cs, xs, 1e-6);                        % 11-by-7
    A  = J(isActive, :).';                              % active constraint gradients
    for i = find(atLower).', e = zeros(7, 1); e(i) = -1; A = [A e]; end %#ok<AGROW>
    for i = find(atUpper).', e = zeros(7, 1); e(i) =  1; A = [A e]; end %#ok<AGROW>
    nActive = size(A, 2);
    if nActive >= 7
        fprintf('%d active constraints and bounds for 7 variables: a vertex, reduced Hessian is empty.\n', nActive);
    else
        Z = null(A.');
        fprintf('Eigenvalues of the reduced Hessian: %s\n', mat2str(sort(eig(Z.'*H*Z)).', 4));
    end

    fid = fopen(fullfile(resDir, ['constraints_' strrep(methods{k}, '-', '_') '.csv']), 'w');
    fprintf(fid, 'constraint,value,limit,unit,margin_pct,active,multiplier\n');
    for i = 1:11
        fprintf(fid, '%s,%.6g,%.6g,%s,%.4f,%d,%.6g\n', conNames{i}, val(i), lim(i), unit{i}, margin(i), isActive(i), b.lamC(i));
    end
    fclose(fid);
end

% Convexity check by sampling: along random line segments inside the bounds,
% a convex function has f(midpoint) <= average of f at the two ends.
fprintf('\n--- Convexity test (midpoint test on 3000 random segments) ---\n');
nSeg = 3000;
rng(2);
nViol = zeros(12, 1);
for s = 1:nSeg
    xa = p.lb + rand(7, 1).*(p.ub - p.lb);
    xb = p.lb + rand(7, 1).*(p.ub - p.lb);
    xm = (xa + xb)/2;
    fa = [objective(xa, p); constraints(xa, p)];
    fb = [objective(xb, p); constraints(xb, p)];
    fm = [objective(xm, p); constraints(xm, p)];
    scale = max(1e-12, (abs(fa) + abs(fb))/2);
    nViol = nViol + ((fm - (fa + fb)/2) > 1e-9*scale);
end
fnNames = [{'objective (volume)'}; conNames];
fid = fopen(fullfile(resDir, 'convexity.csv'), 'w');
fprintf(fid, 'function,violations,segments\n');
for i = 1:12
    fprintf('%-26s violates convexity on %5d of %d segments (%.1f %%)\n', fnNames{i}, nViol(i), nSeg, 100*nViol(i)/nSeg);
    fprintf(fid, '%s,%d,%d\n', fnNames{i}, nViol(i), nSeg);
end
fclose(fid);
fprintf('Sampling can show a function is NOT convex, but it cannot prove convexity.\n\n');

%% 6. REFERENCE CASE (u = 3, T = 1000 Nm) AND GOLINSKI (1970)
% Same limits and bounds as the course case, but the torque and ratio used
% in the validation data. Closer in spirit to Golinski's setting.

fprintf('=== 6. REFERENCE CASE AND LITERATURE ===\n');
pr = p;  pr.T = 1000;  pr.u = 3;
objR = @(x) objective(x, pr);
conR = @(x) constraints(x, pr);
rng(seed);
bestRef = cell(1, 2);
for k = 1:2
    refRuns = [];
    for i = 1:nStartsRef
        x0r = pr.lb + rand(7, 1).*(pr.ub - pr.lb);
        r = run_solver(methods{k}, objR, conR, x0r, pr.lb, pr.ub, 1000);
        if isempty(refRuns), refRuns = r; else, refRuns(end+1) = r; end %#ok<SAGROW>
    end
    bestRef{k} = best_run(refRuns);
    fprintf('%-15s reference case: V = %.4f cm^3, x = %s, success %d/%d\n', methods{k}, ...
        bestRef{k}.f, mat2str(bestRef{k}.x.', 5), sum([refRuns.ok]), nStartsRef);
end

% Golinski Table 1 designs (module converted from cm to mm), evaluated with
% this model at u = 3, T = 1000 Nm
gol = [4.4 6.0 17 7.3 8.1 3.4 5.0;  3.6 7.0 18 6.6 8.2 2.8 5.2].';
golF = [2236.35; 2247.79];
golName = {'Golinski crude Monte Carlo', 'Golinski stray process'};
fprintf('\n%-28s %12s %12s   %s\n', 'design', 'f reported', 'f this model', 'comment');
for i = 1:2
    r = analysis(gol(:, i), pr);
    cv = constraints(gol(:, i), pr);
    outside = varNames(gol(:, i) < pr.lb - 1e-9 | gol(:, i) > pr.ub + 1e-9);
    fprintf('%-28s %12.2f %12.2f   outside our bounds: %s, violated constraints: %d\n', ...
        golName{i}, golF(i), r.V/pr.cm^3, strjoin(outside, ' '), sum(cv > 1e-6));
end
fprintf('%-28s %12s %12.2f   (u = 3, T = 1000 Nm, course limits)\n', 'This work, reference case', '-', bestRef{1}.f);
fprintf('%-28s %12s %12.2f   (u = %g, T = %g Nm)\n', 'This work, course case', '-', best{1}.f, p.u, p.T);
fprintf('Golinski used different loads, limits and bounds, so compare trends only\n');
fprintf('(small z, small module, which constraints are active, variables at bounds).\n\n');

%% 7. FIGURES
if makeFigures
    blue = [0.12 0.37 0.66];  orange = [0.85 0.45 0.10];
    fbest = best{1}.f;

    % convergence history: rerun from start 1 with 1, 2, 3... iterations allowed
    hst = cell(1, 2);
    for k = 1:2
        fullRun = run_solver(methods{k}, obj, con, starts(:, 1), p.lb, p.ub, 1000);
        nIt = max(2, fullRun.iter);
        hst{k} = zeros(nIt, 2);
        for it = 1:nIt
            r = run_solver(methods{k}, obj, con, starts(:, 1), p.lb, p.ub, it);
            hst{k}(it, :) = [abs(r.f - fbest), max(r.viol, 1e-16)];
        end
    end

    % Figure 1: convergence
    f1 = figure('Position', [100 100 800 340]);
    subplot(1, 2, 1); hold on; box on; grid on;
    semilogy(1:size(hst{1}, 1), max(hst{1}(:,1), 1e-9), '-o', 'Color', blue, 'LineWidth', 1.4, 'MarkerSize', 4);
    semilogy(1:size(hst{2}, 1), max(hst{2}(:,1), 1e-9), '-s', 'Color', orange, 'LineWidth', 1.4, 'MarkerSize', 4);
    xlabel('Iteration'); ylabel('|Volume - best known| (cm^3)'); title('(a) Objective error');
    legend('SQP', 'Interior point', 'Location', 'northeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    semilogy(1:size(hst{1}, 1), hst{1}(:,2), '-o', 'Color', blue, 'LineWidth', 1.4, 'MarkerSize', 4);
    semilogy(1:size(hst{2}, 1), hst{2}(:,2), '-s', 'Color', orange, 'LineWidth', 1.4, 'MarkerSize', 4);
    xlabel('Iteration'); ylabel('Max normalised violation (-)'); title('(b) Constraint violation');
    save_fig(f1, figDir, 'fig1_convergence');

    % Figure 2: multistart results
    f2 = figure('Position', [100 100 800 340]);
    subplot(1, 2, 1); hold on; box on; grid on;
    fS = [runs{1}.f] - fbest;  fI = [runs{2}.f] - fbest;
    plot(1:nStarts, fS, 'o', 'Color', blue, 'MarkerSize', 4);
    plot(1:nStarts, fI, 's', 'Color', orange, 'MarkerSize', 4);
    xlabel('Start number'); ylabel('Final volume above best (cm^3)'); title('(a) Final objective per start');
    legend('SQP', 'Interior point', 'Location', 'northeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    edges = linspace(min([fS fI]), max([fS fI]) + 1e-12, 21);
    cS = histc(fS, edges);  cI = histc(fI, edges);
    bar(edges, [cS(:) cI(:)], 'grouped');
    xlabel('Final volume above best (cm^3)'); ylabel('Number of starts'); title('(b) Spread of final values');
    save_fig(f2, figDir, 'fig2_multistart');

    % Figure 3: volume against z
    f3 = figure('Position', [100 100 500 350]); hold on; box on; grid on;
    plot(zList, Vz(:, 1), '-o', 'Color', blue, 'LineWidth', 1.4, 'MarkerSize', 5);
    plot(zList, Vz(:, 2), '--s', 'Color', orange, 'LineWidth', 1.2, 'MarkerSize', 4);
    xlabel('Number of pinion teeth, z (-)'); ylabel('Optimal volume (cm^3)');
    set(gca, 'XTick', zList);
    legend('SQP', 'Interior point', 'Location', 'northwest');
    save_fig(f3, figDir, 'fig3_volume_vs_z');

    % Figure 4: constraint margins at the optimum
    f4 = figure('Position', [100 100 640 460]); hold on; box on; grid on;
    barh(1:11, margins, 'grouped');
    set(gca, 'YTick', 1:11, 'YTickLabel', conNames, 'YDir', 'reverse');
    xlabel('Constraint margin (% of limit), 0 = active');
    legend('SQP', 'Interior point', 'Location', 'southeast');
    save_fig(f4, figDir, 'fig4_constraint_margins');
end

%% SAVE EVERYTHING AND FINISH
save(fullfile(resDir, 'results.mat'), 'p', 'best', 'bestRef', 'zList', 'Vz', 'valComp', 'valErr', 'starts');
fprintf('Done. Results are in %s and figures in %s\n', resDir, figDir);
diary off;


%% LOCAL FUNCTIONS (used by the code above)

function r = analysis(x, p)
% Volume, stresses and deflections for design x (units listed at the top).
x = x(:);
if numel(x) ~= 7 || any(~isfinite(x)) || any(x <= 0)
    error('analysis: x must be 7 positive finite numbers [b m z l1 l2 d1 d2].');
end
% convert to SI
b = x(1)*p.cm;  m = x(2)*p.mm;  z = x(3);
l = [x(4); x(5)]*p.cm;
d = [x(6); x(7)]*p.cm;
T = p.T;  u = p.u;

% gear stresses
r.sigb = p.q*2*T/(b*m^2*z);
r.sigc = sqrt(p.cp*p.kv*T/(b*m^2*z^2)*(1 + u)/u);

% shaft deflection and combined stress (i = 1, 2)
r.y = 8/(3*pi)*T*l.^3./(p.E*d.^4*m*z);
sigsb = 16/pi*T*l./(d.^3*m*z);
sigst = 16/pi*[T; u*T]./d.^3;
r.sigs = sqrt(sigsb.^2 + 3*sigst.^2);

% volume of the two gears and two shafts
zi  = [z; u*z];                          % teeth on gear 1 and gear 2
dst = m*(zi - p.dst_off);
dw  = dst - p.dw_off*m;
dp  = p.dp_fac*d;
Vg = pi*b/4*((dst.^2 - dw.^2) + p.web*(dw.^2 - dp.^2)) + pi/4*p.hub_len*d.*(dp.^2 - d.^2);
Vs = pi/4*d.^2.*l;
r.V = sum(Vg) + sum(Vs);

% other quantities used by the constraints
r.b = b;  r.m = m;  r.l = l;  r.d = d;
r.bm = b/m;
r.size = m*z*(1 + u);
end

function f = objective(x, p)
% Objective: volume divided by 1000 cm^3 (keeps it near 1 for the solver).
r = analysis(x, p);
f = r.V/(1000*p.cm^3);
end

function [c, ceq] = constraints(x, p)
% 11 inequality constraints written as c <= 0.
% Stress, deflection and size rows are value/limit - 1. The four geometric
% rows are linear in x, so they are divided by a constant (the smallest
% allowed length or width) instead of a variable, to keep them linear.
r = analysis(x, p);
refL = p.lb(4)*p.cm;
refB = p.lb(1)*p.cm;
c = [r.sigb/p.sigb_max - 1;
     r.sigc/p.sigc_max - 1;
     r.sigs(1)/p.sigs_max - 1;
     r.sigs(2)/p.sigs_max - 1;
     r.y(1)/p.y_max - 1;
     r.y(2)/p.y_max - 1;
     (p.l1_k*r.d(1) + p.l1_c*p.cm - r.l(1))/refL;
     (p.l2_k*r.d(2) + p.l2_c*p.cm - r.l(2))/refL;
     (p.bm_min*r.m - r.b)/refB;
     (r.b - p.bm_max*r.m)/refB;
     r.size/(p.size_max*p.cm) - 1];
ceq = [];
end

function [val, lim, unit] = physical_constraints(x, p)
% The same 11 constraints in their original physical form (value <= limit),
% used for the tables and for the sanity check.
r = analysis(x, p);
val = [r.sigb/1e6; r.sigc/1e6; r.sigs(1)/1e6; r.sigs(2)/1e6; r.y(1)/p.mm; r.y(2)/p.mm; ...
       p.l1_k*x(6) + p.l1_c; p.l2_k*x(7) + p.l2_c; p.bm_min; r.bm; r.size/p.cm];
lim = [p.sigb_max/1e6; p.sigc_max/1e6; p.sigs_max/1e6; p.sigs_max/1e6; p.y_max/p.mm; p.y_max/p.mm; ...
       x(4); x(5); r.bm; p.bm_max; p.size_max];
unit = {'MPa'; 'MPa'; 'MPa'; 'MPa'; 'mm'; 'mm'; 'cm'; 'cm'; '-'; '-'; 'cm'};
end

function x = insert_z(v, z)
% Build the full design vector from the six free variables and a fixed z.
x = [v(1); v(2); z; v(3); v(4); v(5); v(6)];
end

function r = run_solver(method, obj, con, x0, lb, ub, maxIter)
% One fmincon run with tight tolerances and central finite differences.
opts = optimoptions('fmincon', 'Algorithm', method, 'Display', 'off', ...
    'FiniteDifferenceType', 'central', ...
    'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8, 'StepTolerance', 1e-10, ...
    'MaxIterations', maxIter, 'MaxFunctionEvaluations', 20000);
tic;
[x, fval, flag, output, lambda] = fmincon(obj, x0, [], [], [], [], lb, ub, @(v) nonlcon(con, v), opts);
r.time = toc;
r.x = x;
r.f = fval*1000;                         % volume in cm^3
r.flag = flag;
r.iter = output.iterations;
r.nfev = output.funcCount;
r.firstorder = output.firstorderopt;
r.viol = max([0; con(x)]);               % worst constraint violation
r.ok = (flag > 0) && (r.viol <= 1e-6);   % converged and feasible
r.lamC = lambda.ineqnonlin(:);           % multipliers on the 11 constraints
r.lamLo = lambda.lower(:);               % multipliers on lower bounds
r.lamUp = lambda.upper(:);               % multipliers on upper bounds
end

function [c, ceq] = nonlcon(con, v)
% fmincon wants [c, ceq]
c = con(v);
ceq = [];
end

function b = best_run(runs)
% Lowest volume among the runs that converged and are feasible.
ok = [runs.ok];
if ~any(ok), error('best_run: no run converged to a feasible design.'); end
f = [runs.f];
f(~ok) = inf;
[~, i] = min(f);
b = runs(i);
end

function print_stats(name, runs)
% Multistart statistics for one method.
ok = [runs.ok];
f = [runs.f];
fo = f(ok);
fprintf('\n%s\n', name);
fprintf('  success rate (converged and feasible): %d / %d\n', sum(ok), numel(runs));
fprintf('  best objective  %.6f cm^3\n', min(fo));
fprintf('  worst objective %.6f cm^3\n', max(fo));
fprintf('  spread (worst - best) %.3g cm^3, mean %.6f, std %.3g\n', max(fo) - min(fo), mean(fo), std(fo));
fprintf('  runs within 0.01 %% of the best: %d\n', sum(fo <= min(fo)*(1 + 1e-4)));
fprintf('  iterations: mean %.1f, function evaluations: mean %.1f\n', mean([runs.iter]), mean([runs.nfev]));
fprintf('  run time: mean %.3f s, total %.1f s\n', mean([runs.time]), sum([runs.time]));
end

function save_runs(file, runs)
% Write every multistart run to a csv file.
fid = fopen(file, 'w');
fprintf(fid, 'run,ok,exitflag,V_cm3,max_violation,iterations,func_evals,time_s,b,m,z,l1,l2,d1,d2\n');
for i = 1:numel(runs)
    r = runs(i);
    fprintf(fid, '%d,%d,%d,%.8f,%.3g,%d,%d,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
        i, r.ok, r.flag, r.f, r.viol, r.iter, r.nfev, r.time, r.x);
end
fclose(fid);
end

function g = num_grad(fun, x, h)
% Central difference Jacobian of a scalar or vector function (m-by-n).
n = numel(x);
g = zeros(numel(fun(x)), n);
for j = 1:n
    e = zeros(n, 1);  e(j) = h;
    g(:, j) = (fun(x + e) - fun(x - e))/(2*h);
end
end

function H = num_hessian(fun, x, h)
% Central difference Hessian of a scalar function.
n = numel(x);
H = zeros(n);
for i = 1:n
    for j = 1:n
        ei = zeros(n, 1);  ei(i) = h;
        ej = zeros(n, 1);  ej(j) = h;
        H(i, j) = (fun(x + ei + ej) - fun(x + ei - ej) - fun(x - ei + ej) + fun(x - ei - ej))/(4*h^2);
    end
end
H = (H + H.')/2;
end

function save_fig(f, folder, name)
% Save a figure as PNG (300 dpi) and PDF.
set(f, 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
if exist('exportgraphics', 'file') ~= 0
    exportgraphics(f, fullfile(folder, [name '.png']), 'Resolution', 300);
    exportgraphics(f, fullfile(folder, [name '.pdf']));
else
    print(f, fullfile(folder, [name '.png']), '-dpng', '-r300');
    print(f, fullfile(folder, [name '.pdf']), '-dpdf');
end
close(f);
end
