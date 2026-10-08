% VALIDATE_SPEED_REDUCER  Check the analysis against the brief's validation data.
%
% Author: AWD Labs
% Student ID: 52104479
%
% Runs T = 1000 Nm and T = 2400 Nm, prints computed vs target with percent
% errors, and tests the two readings of the elastic coefficient cp.
% Writes results/matlab/validation.csv. Run from code/matlab.

clear; clc;
xVal = [3.5; 7.0; 22; 7.4; 7.8; 3.5; 5.2];
names = {'Volume (cm^3)','sigma_b (MPa)','sigma_c (MPa)','sigma_s1 (MPa)', ...
         'sigma_s2 (MPa)','y1 (mm)','y2 (mm)'};
target = [4147 323 532 512 454 0.0179 0.0043];

cases = {'validation1000','validation2400'};
comp = zeros(numel(target), 2);
for k = 1:2
    a = speed_reducer_analysis(xVal, speed_reducer_params(cases{k}));
    d = a.display;
    comp(:,k) = [d.V_cm3; d.sigma_b_MPa; d.sigma_c_MPa; d.sigma_s_MPa(1); ...
                 d.sigma_s_MPa(2); d.y_mm(1); d.y_mm(2)];
end
err = 100*(comp - target(:))./target(:);

fprintf('%-16s %10s | %10s %8s | %10s %8s\n', 'quantity','target','T=1000','err %','T=2400','err %');
for i = 1:numel(target)
    fprintf('%-16s %10.4g | %10.4g %8.2f | %10.4g %8.2f\n', names{i}, target(i), ...
            comp(i,1), err(i,1), comp(i,2), err(i,2));
end

fprintf('\nsigma_c readings at T = 2400 Nm (target 532 MPa):\n');
forms = {'one_minus_nu2','one_plus_nu2'};
for k = 1:2
    pk = speed_reducer_params('validation2400', struct('cp_form', forms{k}));
    sc = speed_reducer_analysis(xVal, pk).display.sigma_c_MPa;
    fprintf('  %-14s sigma_c = %.2f MPa (err %.2f %%)\n', forms{k}, sc, 100*(sc-532)/532);
end

outDir = fullfile('..','..','results','matlab');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fid = fopen(fullfile(outDir,'validation.csv'), 'w');
fprintf(fid, 'quantity,target,T1000_computed,T1000_err_pct,T2400_computed,T2400_err_pct\n');
for i = 1:numel(target)
    fprintf(fid, '%s,%g,%g,%g,%g,%g\n', names{i}, target(i), comp(i,1), err(i,1), comp(i,2), err(i,2));
end
fclose(fid);
