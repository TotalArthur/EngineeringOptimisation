% RUN_OPTIMISATION  Main driver: multistart, integer teeth study, reference case.
%
% Author: AWD Labs
% Student ID: 52104479
%
% Run from code/matlab. Requires the Optimization Toolbox.
% Saves results/matlab/*.mat and *.csv. Plot and post processing scripts are
% separate (postprocess_design.m, make_figures.m).
% Random starts: rng(1), uniform in the scaled box [0,1]^n.

clear; clc;
methods = {'sqp', 'interior-point'};
nStarts = 100;     % per method, continuous multistart
nStartsZ = 20;     % per method and per z
nStartsRef = 50;
outDir = fullfile('..','..','results','matlab');
if ~exist(outDir, 'dir'), mkdir(outDir); end

validate_speed_reducer;   % step 1: validation (prints table)

p = speed_reducer_params('optimisation');
fprintf('\nOptimisation case: T = %g Nm, u = %g (student ID %s)\n', p.T, p.u, p.student_id);

%% 2. continuous multistart
rng(1);
starts = rand(7, nStarts);
ms = struct();
for im = 1:2
    m = methods{im};
    runs = cell(nStarts, 1);
    for k = 1:nStarts
        runs{k} = solve_speed_reducer(m, starts(:,k), p, [], false);
    end
    ms.(matlab.lang.makeValidName(m)) = runs;
    printStats(m, runs);
    writeRunsCsv(fullfile(outDir, ['multistart_' strrep(m,'-','_') '.csv']), runs);
end
best = struct();
for im = 1:2
    m = matlab.lang.makeValidName(methods{im});
    best.(m) = bestOf(ms.(m));
end
% convergence history of the first start for both methods
hist1 = struct();
for im = 1:2
    r = solve_speed_reducer(methods{im}, starts(:,1), p, [], true);
    hist1.(matlab.lang.makeValidName(methods{im})) = r.history;
end
save(fullfile(outDir, 'multistart.mat'), 'starts', 'ms', 'best', 'hist1', 'p');

%% 3. integer teeth study
zs = 17:28;
Vz = nan(numel(zs), 2);
Xz = nan(7, numel(zs), 2);
for iz = 1:numel(zs)
    fixed = struct('idx', p.z_index, 'val', zs(iz));
    rng(1);
    startsZ = rand(6, nStartsZ);
    for im = 1:2
        runs = cell(nStartsZ, 1);
        for k = 1:nStartsZ
            runs{k} = solve_speed_reducer(methods{im}, startsZ(:,k), p, fixed, false);
        end
        b = bestOf(runs);
        if ~isempty(b), Vz(iz,im) = b.f; Xz(:,iz,im) = b.x; end
    end
    fprintf('z = %d: SQP %.4f, IP %.4f cm^3\n', zs(iz), Vz(iz,1), Vz(iz,2));
end
save(fullfile(outDir, 'integer_z.mat'), 'zs', 'Vz', 'Xz');
writematrix([zs(:), Vz], fullfile(outDir, 'integer_z.csv'));

%% 4. reference case u = 3, T = 1000 Nm
pr = speed_reducer_params('reference');
rng(1);
startsR = rand(7, nStartsRef);
bestRef = struct();
for im = 1:2
    runs = cell(nStartsRef, 1);
    for k = 1:nStartsRef
        runs{k} = solve_speed_reducer(methods{im}, startsR(:,k), pr, [], false);
    end
    bestRef.(matlab.lang.makeValidName(methods{im})) = bestOf(runs);
    printStats(['reference ' methods{im}], runs);
end
save(fullfile(outDir, 'reference.mat'), 'bestRef', 'pr');

%% local functions
function b = bestOf(runs)
ok = cellfun(@(r) r.ok, runs);
b = [];
if any(ok)
    f = cellfun(@(r) r.f, runs);
    f(~ok) = inf;
    [~, i] = min(f);
    b = runs{i};
end
end

function printStats(name, runs)
ok = cellfun(@(r) r.ok, runs);
f = cellfun(@(r) r.f, runs);
fprintf('%-22s success %d/%d', name, sum(ok), numel(runs));
if any(ok)
    fo = f(ok);
    fprintf(' best %.6f worst %.6f spread %.3g mean it %.1f mean nfev %.1f mean time %.3f s\n', ...
        min(fo), max(fo), max(fo)-min(fo), mean(cellfun(@(r) r.nit, runs)), ...
        mean(cellfun(@(r) r.nfev, runs)), mean(cellfun(@(r) r.time, runs)));
else
    fprintf('\n');
end
end

function writeRunsCsv(path, runs)
fid = fopen(path, 'w');
fprintf(fid, 'run,ok,exitflag,f_cm3,max_viol,iterations,func_evals,time_s,b_cm,m_mm,z,l1_cm,l2_cm,d1_cm,d2_cm\n');
for k = 1:numel(runs)
    r = runs{k};
    fprintf(fid, '%d,%d,%d,%.10g,%.3g,%d,%d,%.4f,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n', ...
        k, r.ok, r.exitflag, r.f, r.viol, r.nit, r.nfev, r.time, r.x);
end
fclose(fid);
end
