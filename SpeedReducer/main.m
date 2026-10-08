% main.m
% Author: AWD Labs
% EG503X Engineering Optimisation, Assignment 1
% Speed reducer: minimise the material volume of two gears and two shafts
%
% HOW TO RUN
%   1. Keep all six .m files in the same folder (main, speedReducerAnalysis,
%      objectiveFun, constraintFun, solveDesign, printConstraintTable).
%   2. Check the student ID digits in Section 1, then press Run.
%   3. Results print to the command window, figures are saved as PNG and a
%      copy of the printout is saved, all in a folder called "results".
%
% Design vector used everywhere (brief units, so it is well scaled):
%   x = [b (cm), m (mm), z (-), l1 (cm), l2 (cm), d1 (cm), d2 (cm)]

clear; clc; close all;

%% 1. STUDENT ID, TOOLBOX CHECK AND ALL PARAMETERS
% Student ID 52104479, so the last digit is 9 and the second last is 7.
n    = 9;      % last digit of student ID
m_id = 7;      % second last digit of student ID

if isnan(n) || isnan(m_id)
    error('main:noID', ['Student ID digits are missing. Set n (last digit) ' ...
          'and m_id (second last digit) at the top of main.m, then run again.']);
end

% fmincon lives in the Optimization Toolbox, so stop early if it is missing.
if ~license('test', 'Optimization_Toolbox') || isempty(which('fmincon'))
    error('main:noToolbox', ['The Optimization Toolbox is not installed or not ' ...
          'licensed, so fmincon is not available. Install it from Home > Add-Ons ' ...
          'and run again.']);
end

% Output folder for figures and the saved printout
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
outDir = fullfile(thisDir, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
diary(fullfile(outDir, 'command_window_output.txt'));   % copy of everything printed

% ---- Everything below is in SI units (m, Nm, Pa, m^3) unless stated ----
params.n    = n;
params.m_id = m_id;
params.T    = 1500 + 100*n;        % input torque (Nm)
params.u    = 2.5 + 0.1*m_id;      % speed ratio (-)
params.q    = 2.54;                % gear tooth form factor (-)
params.kv   = 2.1;                 % contact stress concentration factor (-)
params.E    = 200e9;               % Young's modulus (Pa)
params.nu   = 0.3;                 % Poisson's ratio (-)
params.cpNuSign = -1;              % -1 gives cp = E/(2 pi (1 - nu^2)), +1 gives (1 + nu^2)

% Mechanical limits
params.sigma_b_max = 650e6;        % gear bending stress (Pa)
params.sigma_c_max = 800e6;        % gear contact stress (Pa)
params.sigma_s_max = 550e6;        % combined shaft stress (Pa)
params.y_max       = 0.075e-3;     % shaft deflection (m), brief says 0.075 mm

% Geometric requirements
params.l1_perD1 = 1.5;             % l1 >= 1.5 d1 + 1.9 cm
params.l1_extra = 1.9e-2;          % 1.9 cm in m
params.l2_perD2 = 1.1;             % l2 >= 1.1 d2 + 1.9 cm
params.l2_extra = 1.9e-2;
params.bm_min   = 5;               % 5 <= b/m <= 12
params.bm_max   = 12;
params.size_max = 160e-2;          % m z (1+u) <= 160 cm, in m

% Gear geometry constants used in the volume equations (Figure 2 of the brief)
params.kTip      = 2.4;            % d_st = m (z - 2.4)
params.kRim      = 4;              % d_w  = d_st - 4 m
params.kHub      = 2.4;            % d_p  = 2.4 d
params.kWeb      = 1/3;            % web width is b/3
params.kHubWidth = 2;              % hub width is 2 d

% Design variable bounds, in the brief's units [b m z l1 l2 d1 d2]
params.lb = [2.6, 7.0, 17, 7.3, 7.3, 2.8, 5.0];
params.ub = [4.4, 8.0, 28, 8.3, 8.3, 3.9, 5.5];
params.xToSI   = [1e-2, 1e-3, 1, 1e-2, 1e-2, 1e-2, 1e-2];  % brief units -> SI
params.xNames  = {'b', 'm', 'z', 'l1', 'l2', 'd1', 'd2'};
params.xUnits  = {'cm', 'mm', '-', 'cm', 'cm', 'cm', 'cm'};

% Constraint list (same order as r.conValue in speedReducerAnalysis)
params.conNames  = {'Gear bending stress', 'Gear contact stress', ...
                    'Shaft 1 combined stress', 'Shaft 2 combined stress', ...
                    'Shaft 1 deflection', 'Shaft 2 deflection', ...
                    'Shaft 1 length (min)', 'Shaft 2 length (min)', ...
                    'Face width / module (min)', 'Face width / module (max)', ...
                    'Overall size m z (1+u)'};
params.conUnits  = {'MPa', 'MPa', 'MPa', 'MPa', 'mm', 'mm', 'cm', 'cm', '-', '-', 'cm'};
params.conToDisp = [1e-6, 1e-6, 1e-6, 1e-6, 1e3, 1e3, 1e2, 1e2, 1, 1, 1e2]; % SI -> printed units
params.conIsMin  = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0];   % 1 = value must stay ABOVE the limit

% Unit factors for printing results in the brief's units
params.toMPa = 1e-6;               % Pa  -> MPa
params.toMM  = 1e3;                % m   -> mm
params.toCm3 = 1e6;                % m^3 -> cm^3
params.objScale = 1000;            % objective is volume (cm^3) / 1000

% Solver settings (same for both algorithms so the comparison is fair)
params.optTol  = 1e-7;
params.conTol  = 1e-9;
params.stepTol = 1e-10;
params.maxIter = 1000;
params.maxEval = 20000;

% Tolerances used when reporting
params.feasTol   = 1e-6;           % c <= feasTol counts as feasible (0.0001 % of the limit)
params.activeTol = 1e-3;           % utilisation >= 99.9 % counts as ACTIVE
params.boundTol  = 1e-4;           % within 0.01 % of the bound range counts as ON the bound
params.sameTol   = 1e-3;           % volumes within 0.1 % count as the same optimum
params.nStarts   = 20;             % number of random starting points
params.seed      = 503;            % random seed so the runs repeat

algs  = {'sqp', 'interior-point'};
nAlg  = numel(algs);
nCon  = numel(params.conNames);

fprintf('==================================================================\n');
fprintf(' SPEED REDUCER OPTIMISATION   (EG503X Assignment 1)\n');
fprintf('==================================================================\n');
fprintf('Student ID digits: n = %d, m = %d\n', n, m_id);
fprintf('Input torque T = 1500 + 100*%d = %g Nm\n', n, params.T);
fprintf('Speed ratio  u = 2.5 + 0.1*%d = %.2f\n\n', m_id, params.u);

%% 2. VALIDATION CASE (brief, page 4)
fprintf('------------------------------------------------------------------\n');
fprintf(' 2. VALIDATION CASE\n');
fprintf('------------------------------------------------------------------\n');

pv = params;                 % copy of params with the validation inputs
pv.u = 3;
pv.T = 1000;                 % torque as written in the brief
xVal = [3.5, 7.0, 22, 7.4, 7.8, 3.5, 5.2];   % [b m z l1 l2 d1 d2] in brief units
expected = [4147, 323, 532, 512, 454, 0.0179, 0.0043];   % from the brief
valNames = {'Volume (cm^3)', 'sigma_b (MPa)', 'sigma_c (MPa)', ...
            'sigma_s1 (MPa)', 'sigma_s2 (MPa)', 'y1 (mm)', 'y2 (mm)'};

% anonymous function that turns the analysis struct into the 7 numbers in the table
getNumbers = @(r) [r.V_total*params.toCm3, r.sigma_b*params.toMPa, ...
                   r.sigma_c*params.toMPa, r.sigma_s(1)*params.toMPa, ...
                   r.sigma_s(2)*params.toMPa, r.y(1)*params.toMM, r.y(2)*params.toMM];

calcAsStated = getNumbers(speedReducerAnalysis(xVal, pv));

% The stresses and deflections all come out about 2.4 times too small. Every
% one of them is proportional to T (except sigma_c, which goes with sqrt(T)),
% so the torque is the suspect. Work out which torque would give the expected
% bending stress, then check whether that ONE number fixes everything else.
T_implied = pv.T * expected(2) / calcAsStated(2);
T_diag    = round(T_implied/100) * 100;
pd = pv;
pd.T = T_diag;
calcDiag = getNumbers(speedReducerAnalysis(xVal, pd));

diffAsStated = (calcAsStated - expected) ./ expected * 100;
diffDiag     = (calcDiag     - expected) ./ expected * 100;

fprintf('Inputs: u = %g, z = %g, m = %g mm, b = %g cm, l1 = %g, l2 = %g, d1 = %g, d2 = %g cm\n\n', ...
        pv.u, xVal(3), xVal(2), xVal(1), xVal(4), xVal(5), xVal(6), xVal(7));
fprintf('%-16s %10s | %12s %9s %9s | %12s %9s\n', 'Quantity', 'Expected', ...
        sprintf('Calc T=%g', pv.T), '% diff', 'Exp/Calc', sprintf('Calc T=%g', T_diag), '% diff');
for k = 1:7
    fprintf('%-16s %10.4f | %12.4f %9.2f %9.3f | %12.4f %9.2f\n', valNames{k}, ...
            expected(k), calcAsStated(k), diffAsStated(k), expected(k)/calcAsStated(k), ...
            calcDiag(k), diffDiag(k));
end

fprintf('\nTorque implied by the expected bending stress: %.1f Nm (rounded to %g Nm)\n', T_implied, T_diag);
fprintf('Largest difference in stresses and deflections at T = %g Nm: %.2f %%\n', ...
        T_diag, max(abs(diffDiag(2:7))));
fprintf('(expected values in the brief are rounded to 3 significant figures)\n');

% Which version of the contact stress equation matches?
fprintf('\nContact stress sigma_c (MPa), expected %g:\n', expected(3));
fprintf('%-52s %10s %10s\n', 'Version of equation', sprintf('T=%g', pv.T), sprintf('T=%g', T_diag));
signNames = {'cp = E/(2 pi (1 + nu^2)), as printed in the brief', ...
             'cp = E/(2 pi (1 - nu^2)), standard Hertz contact'};
signValues = [1, -1];
for s = 1:2
    row = zeros(1, 2);
    torques = [pv.T, T_diag];
    for t = 1:2
        pt = pv;
        pt.T = torques(t);
        pt.cpNuSign = signValues(s);
        rTmp = speedReducerAnalysis(xVal, pt);
        row(t) = rTmp.sigma_c * params.toMPa;
    end
    fprintf('%-52s %10.1f %10.1f\n', signNames{s}, row(1), row(2));
end
% The equation as it was typed in my task notes: kv T/(b m^2 z) * 2 (1+u)/u, no z^2
xv = xVal .* params.xToSI;
cpMinus = pv.E / (2*pi*(1 - pv.nu^2));
sigmaAlt = sqrt(cpMinus * pv.kv * T_diag / (xv(1)*xv(2)^2*xv(3)) * 2*(1+pv.u)/pv.u);
fprintf('%-52s %10s %10.1f\n', 'Notes version (no z^2, factor 2), (1 - nu^2)', '-', sigmaAlt * params.toMPa);

% Volume check
rDiag = speedReducerAnalysis(xVal, pd);
fprintf('\nVolume check (does not depend on T):\n');
fprintf('  gears  = %.1f cm^3 (gear 1 %.1f, gear 2 %.1f)\n', rDiag.V_gear_total*params.toCm3, ...
        rDiag.V_gears(1)*params.toCm3, rDiag.V_gears(2)*params.toCm3);
fprintf('  shafts = %.1f cm^3 (cylinders, pi/4 d^2 l)\n', rDiag.V_shaft_total*params.toCm3);
fprintf('  total  = %.1f cm^3 against %.0f cm^3 expected, difference %.1f cm^3 (%.2f %%)\n', ...
        rDiag.V_total*params.toCm3, expected(1), rDiag.V_total*params.toCm3 - expected(1), diffDiag(1));
fprintf('  The volume equations as written in the brief do not reproduce 4147 cm^3.\n');
fprintf('  See the notes in the report. The optimum below uses the equations as given.\n\n');

%% 3. STARTING POINT AND CONSTRAINTS AT THE START
fprintf('------------------------------------------------------------------\n');
fprintf(' 3. STARTING POINT (feasible)\n');
fprintf('------------------------------------------------------------------\n');
x0 = [4.0, 7.5, 22, 7.8, 7.8, 3.6, 5.2];   % [b m z l1 l2 d1 d2]
for k = 1:7
    fprintf('  %-3s = %7.3f %s\n', params.xNames{k}, x0(k), params.xUnits{k});
end
fprintf('  Volume at start = %.1f cm^3\n\n', objectiveFun(x0, params) * params.objScale);
printConstraintTable(x0, params);
if max(constraintFun(x0, params)) > 0
    warning('The starting point is not feasible. Change x0 in Section 3.');
end
fprintf('\n');

%% 4. OPTIMISATION WITH TWO ALGORITHMS FROM THE SAME START
fprintf('------------------------------------------------------------------\n');
fprintf(' 4. fmincon: SQP versus INTERIOR-POINT (same start, same options)\n');
fprintf('------------------------------------------------------------------\n');

% one throw-away run per algorithm so MATLAB loads everything before the timed runs
for a = 1:nAlg
    solveDesign(algs{a}, x0, params);
end

for a = 1:nAlg
    sol(a) = solveDesign(algs{a}, x0, params);
end

fprintf('%-6s %-5s %16s %16s\n', 'Var', 'Unit', algs{1}, algs{2});
for k = 1:7
    fprintf('%-6s %-5s %16.4f %16.4f\n', params.xNames{k}, params.xUnits{k}, ...
            sol(1).x(k), sol(2).x(k));
end
fprintf('\n%-28s %16.2f %16.2f\n', 'Volume (cm^3)',         sol(1).V, sol(2).V);
fprintf('%-28s %16d %16d\n',       'Exit flag',             sol(1).exitflag, sol(2).exitflag);
fprintf('%-28s %16d %16d\n',       'Iterations',            sol(1).iterations, sol(2).iterations);
fprintf('%-28s %16d %16d\n',       'Function evaluations',  sol(1).funcCount, sol(2).funcCount);
fprintf('%-28s %16.3f %16.3f\n',   'Run time (s)',          sol(1).time, sol(2).time);
fprintf('%-28s %16.2e %16.2e\n',   'Largest constraint c',  sol(1).maxViol, sol(2).maxViol);
fprintf('(exit flag > 0 means fmincon converged; c <= 0 means feasible)\n\n');

%% 5. CONSTRAINT TABLES AND BOUNDS AT THE OPTIMUM
fprintf('------------------------------------------------------------------\n');
fprintf(' 5. CONSTRAINTS AND BOUNDS AT THE OPTIMUM\n');
fprintf('------------------------------------------------------------------\n');
fprintf('ACTIVE = utilisation within %.1g of 100 %%.  Bound test tolerance = %.1g of the bound range.\n\n', ...
        params.activeTol, params.boundTol);

for a = 1:nAlg
    fprintf('--- %s ---\n', algs{a});
    printConstraintTable(sol(a).x, params);
    fprintf('\nDesign variables on a bound (%s):\n', algs{a});
    nOnBound = 0;
    for k = 1:7
        range = params.ub(k) - params.lb(k);
        if abs(sol(a).x(k) - params.lb(k)) <= params.boundTol * range
            bndText = 'ON LOWER BOUND';
            nOnBound = nOnBound + 1;
        elseif abs(sol(a).x(k) - params.ub(k)) <= params.boundTol * range
            bndText = 'ON UPPER BOUND';
            nOnBound = nOnBound + 1;
        else
            bndText = 'free';
        end
        fprintf('  %-3s = %8.4f %-3s  [%5.2f, %5.2f]  %s\n', params.xNames{k}, ...
                sol(a).x(k), params.xUnits{k}, params.lb(k), params.ub(k), bndText);
    end
    fprintf('  %d of 7 variables sit on a bound.\n\n', nOnBound);
end

%% 6. MULTI-START: 20 RANDOM STARTING POINTS FOR EACH METHOD
fprintf('------------------------------------------------------------------\n');
fprintf(' 6. MULTI-START (%d random starts, seed %d)\n', params.nStarts, params.seed);
fprintf('------------------------------------------------------------------\n');

rng(params.seed);
nS = params.nStarts;
lbRows = repmat(params.lb, nS, 1);
ubRows = repmat(params.ub, nS, 1);
X0 = lbRows + rand(nS, 7) .* (ubRows - lbRows);     % uniform inside the bounds

V_ms    = nan(nS, nAlg);      % final volume of each run (cm^3)
flag_ms = nan(nS, nAlg);      % exit flag of each run
good_ms = false(nS, nAlg);    % true if converged AND feasible
for k = 1:nS
    for a = 1:nAlg
        s = solveDesign(algs{a}, X0(k, :), params);
        V_ms(k, a)    = s.V;
        flag_ms(k, a) = s.exitflag;
        good_ms(k, a) = (s.exitflag > 0) && (s.maxViol <= params.feasTol);
    end
end

Vbest = min(V_ms(good_ms));       % best volume found by any good run
if isempty(Vbest), Vbest = NaN; end
fprintf('Best volume found by any run: %.2f cm^3\n\n', Vbest);
fprintf('%-16s %9s %12s %14s %10s %10s %10s %10s\n', 'Method', 'Starts', ...
        'Converged', 'Same optimum', 'Min V', 'Max V', 'Spread', 'Std dev');
for a = 1:nAlg
    Vg = V_ms(good_ms(:, a), a);
    nSame = sum(abs(Vg - Vbest) <= params.sameTol * Vbest);
    if isempty(Vg)
        fprintf('%-16s %9d %12d %14d   (no successful runs)\n', algs{a}, nS, 0, 0);
    else
        fprintf('%-16s %9d %12d %14d %10.2f %10.2f %10.4f %10.4f\n', algs{a}, nS, ...
                numel(Vg), nSame, min(Vg), max(Vg), max(Vg) - min(Vg), std(Vg));
    end
end
fprintf('(Converged = exit flag > 0 and feasible. Same optimum = within %.1g of the best volume.)\n\n', ...
        params.sameTol);

fprintf('Run-by-run volumes (cm^3), exit flags in brackets:\n');
fprintf('%5s %20s %20s\n', 'Run', algs{1}, algs{2});
for k = 1:nS
    fprintf('%5d %14.2f (%2d) %14.2f (%2d)\n', k, V_ms(k, 1), flag_ms(k, 1), ...
            V_ms(k, 2), flag_ms(k, 2));
end
fprintf('\n');

%% 7. ROUND z TO AN INTEGER AND RE-SOLVE THE OTHER VARIABLES
fprintf('------------------------------------------------------------------\n');
fprintf(' 7. ROUNDED z (z fixed, other six variables re-optimised)\n');
fprintf('------------------------------------------------------------------\n');
for a = 1:nAlg
    zCont  = sol(a).x(3);
    zRound = round(zCont);
    solR(a) = solveDesign(algs{a}, sol(a).x, params, zRound);
    dV = solR(a).V - sol(a).V;
    fprintf('--- %s ---\n', algs{a});
    fprintf('  z continuous = %.4f, z rounded = %d (gear 2 would have u*z = %.2f teeth)\n', ...
            zCont, zRound, params.u * zRound);
    fprintf('  Volume: continuous %.2f cm^3, rounded %.2f cm^3, change %+.3f cm^3 (%+.4f %%)\n', ...
            sol(a).V, solR(a).V, dV, dV / sol(a).V * 100);
    if solR(a).maxViol <= params.feasTol, feasText = 'YES'; else feasText = 'NO'; end
    fprintf('  Feasible: %s (largest constraint c = %.2e), exit flag %d\n', ...
            feasText, solR(a).maxViol, solR(a).exitflag);
    fprintf('  x = [');
    for k = 1:7
        fprintf(' %s=%.4f', params.xNames{k}, solR(a).x(k));
    end
    fprintf(' ]\n\n');
end

%% 8. FIGURES (saved as PNG in the results folder)
% Figure 1: how close each constraint is to its limit at the optimum
util = zeros(nCon, nAlg);
for a = 1:nAlg
    util(:, a) = (constraintFun(sol(a).x, params)' + 1) * 100;
end
fig1 = figure('Name', 'Constraint utilisation', 'Position', [100 100 900 550]);
barh(util);
hold on;
plot([100 100], [0.4, nCon + 0.6], 'r--', 'LineWidth', 1.5);
set(gca, 'YTick', 1:nCon, 'YTickLabel', params.conNames, 'YDir', 'reverse');
xlabel('Constraint utilisation (% of limit)');
ylabel('Constraint');
title('Constraint utilisation at the optimum (100 % = on the limit)');
xlim([0 140]);                 % leave room on the right for the legend
legend([algs, {'Limit (100 %)'}], 'Location', 'southeast');
grid on;
print(fig1, fullfile(outDir, 'fig1_constraint_utilisation.png'), '-dpng', '-r200');

% Figure 2: final volume for each random start
fig2 = figure('Name', 'Multi-start volumes', 'Position', [100 100 800 500]);
plot(1:nS, V_ms(:, 1), 'o-', 'LineWidth', 1.2, 'MarkerSize', 7);
hold on;
plot(1:nS, V_ms(:, 2), 's--', 'LineWidth', 1.2, 'MarkerSize', 7);
for a = 1:nAlg
    bad = ~good_ms(:, a);
    if any(bad)
        plot(find(bad), V_ms(bad, a), 'rx', 'MarkerSize', 12, 'LineWidth', 2, ...
             'HandleVisibility', 'off');
    end
end
xlabel('Random start number');
ylabel('Final volume (cm^3)');
title(sprintf('Final volume from %d random starts (red x = not converged or infeasible)', nS));
% +/- 5 % window, otherwise MATLAB zooms in on tiny numerical noise when every
% run lands on the same volume and makes it look like a big spread
ylim([min(V_ms(:)) * 0.95, max(V_ms(:)) * 1.05]);
legend(algs, 'Location', 'northeast');
grid on;
print(fig2, fullfile(outDir, 'fig2_multistart_volume.png'), '-dpng', '-r200');

% Figure 3: where each design variable sits between its bounds
pos = zeros(7, nAlg);
for a = 1:nAlg
    pos(:, a) = ((sol(a).x - params.lb) ./ (params.ub - params.lb) * 100)';
end
varLabels = cell(1, 7);
for k = 1:7
    varLabels{k} = sprintf('%s (%s)', params.xNames{k}, params.xUnits{k});
end
fig3 = figure('Name', 'Variables within bounds', 'Position', [100 100 800 500]);
bar(pos);
set(gca, 'XTick', 1:7, 'XTickLabel', varLabels);
ylim([0 105]);
xlabel('Design variable');
ylabel('Position within bounds (0 % = lower, 100 % = upper)');
title('Design variables at the optimum relative to their bounds');
legend(algs, 'Location', 'northeast');
grid on;
print(fig3, fullfile(outDir, 'fig3_variables_vs_bounds.png'), '-dpng', '-r200');

fprintf('Figures and printout saved in: %s\n', outDir);
diary off;
