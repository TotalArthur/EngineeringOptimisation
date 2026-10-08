% MAKE_FIGURES  Report figures (PNG 300 dpi + PDF) from results/matlab.
%
% Author: AWD Labs
% Student ID: 52104479
%
% Separate from the solver. Blue/orange pair with different markers so the
% two methods never rely on colour alone. Run after run_optimisation and
% postprocess_design.

clear; clc;
outDir = fullfile('..','..','results','matlab');
figDir = fullfile('..','..','figures','matlab');
if ~exist(figDir, 'dir'), mkdir(figDir); end
S = load(fullfile(outDir, 'multistart.mat'));
Z = load(fullfile(outDir, 'integer_z.mat'));
P = load(fullfile(outDir, 'postprocess.mat'));
blue = [31 95 168]/255;  orange = [217 115 26]/255;
col = {blue, orange};  mk = {'o', 's'};  lbl = {'SQP', 'Interior point'};
fn = {'sqp', 'interior_point'};
set(groot, 'defaultAxesFontSize', 11, 'defaultTextFontSize', 11);

fbest = S.best.sqp.f;

% Fig 1 convergence
f = figure('Position', [100 100 800 340]);
subplot(1,2,1); hold on; box on; grid on;
for k = 1:2
    h = S.hist1.(fn{k});
    semilogy(1:size(h,1), max(abs(h(:,1) - fbest), 1e-9), ['-' mk{k}], 'Color', col{k}, 'MarkerSize', 4, 'LineWidth', 1.4);
end
set(gca, 'YScale', 'log'); xlabel('Iteration'); ylabel('|Volume - best known| (cm^3)'); title('(a) Objective error');
legend(lbl, 'Location', 'northeast', 'Box', 'off');
subplot(1,2,2); hold on; box on; grid on;
for k = 1:2
    h = S.hist1.(fn{k});
    semilogy(1:size(h,1), max(h(:,2), 1e-16), ['-' mk{k}], 'Color', col{k}, 'MarkerSize', 4, 'LineWidth', 1.4);
end
set(gca, 'YScale', 'log'); xlabel('Iteration'); ylabel('Max normalised violation (-)'); title('(b) Constraint violation');
saveBoth(f, figDir, 'fig1_convergence');

% Fig 2 multistart
f = figure('Position', [100 100 800 340]);
subplot(1,2,1); hold on; box on; grid on;
fall = cell(1,2);
for k = 1:2
    runs = S.ms.(fn{k});
    fall{k} = cellfun(@(r) r.f, runs);
    ok = cellfun(@(r) r.ok, runs);
    plot(find(ok), fall{k}(ok) - fbest, mk{k}, 'Color', col{k}, 'MarkerSize', 4);
    plot(find(~ok), fall{k}(~ok) - fbest, 'kx', 'MarkerSize', 6);
end
xlabel('Start number'); ylabel('Final volume above best (cm^3)'); title('(a) Final objective per start');
legend(lbl, 'Location', 'northeast', 'Box', 'off');
subplot(1,2,2); hold on; box on; grid on;
for k = 1:2
    histogram(fall{k} - fbest, 20, 'FaceColor', col{k}, 'FaceAlpha', 0.6);
end
xlabel('Final volume above best (cm^3)'); ylabel('Number of starts'); title('(b) Spread of final values');
saveBoth(f, figDir, 'fig2_multistart');

% Fig 3 volume vs z
f = figure('Position', [100 100 500 350]); hold on; box on; grid on;
plot(Z.zs, Z.Vz(:,1), '-o', 'Color', blue, 'LineWidth', 1.4, 'MarkerSize', 5);
plot(Z.zs, Z.Vz(:,2), '--s', 'Color', orange, 'LineWidth', 1.2, 'MarkerSize', 4);
xlabel('Number of pinion teeth, z (-)'); ylabel('Optimal volume (cm^3)'); xticks(Z.zs);
legend(lbl, 'Location', 'northwest', 'Box', 'off');
saveBoth(f, figDir, 'fig3_volume_vs_z');

% Fig 4 constraint margins
f = figure('Position', [100 100 640 440]); hold on; box on; grid on;
names = P.post.sqp.constraints.name;
y = (1:numel(names)).';
barh(y, [P.post.sqp.margin_pct, P.post.interior_point.margin_pct], 'grouped');
colormap([blue; orange]);
set(gca, 'YTick', y, 'YTickLabel', names, 'YDir', 'reverse');
xlabel('Constraint margin (% of limit), 0 = active');
legend(lbl, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
saveBoth(f, figDir, 'fig4_constraint_margins');

function saveBoth(f, d, name)
set(f, 'PaperPositionMode', 'auto');
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(d, [name '.pdf']), 'ContentType', 'vector');
close(f);
end
