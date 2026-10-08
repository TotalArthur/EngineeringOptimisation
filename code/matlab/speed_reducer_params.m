function p = speed_reducer_params(caseName, overrides, studentId)
% SPEED_REDUCER_PARAMS  Every constant, limit and bound for the speed reducer.
%
% Author: AWD Labs
% Student ID: 52104479
% Course: EG503X Engineering Optimisation, Assignment 1
%
% UNIT CONVENTION (whole code base)
%   Design vector x is in DISPLAY units: [b cm, m mm, z, l1 cm, l2 cm, d1 cm, d2 cm].
%   speed_reducer_analysis converts to SI (m, Pa, N m, m^3) immediately and
%   works in SI. Results are converted back only for display.
%
% p = speed_reducer_params(caseName, overrides, studentId)
%   caseName : 'optimisation' (default), 'validation1000', 'validation2400', 'reference'
%   overrides: struct of fields to replace, e.g. struct('cp_form','one_plus_nu2')
%   studentId: char, default '52104479'. Last digit n gives T = 1500 + 100 n,
%              second last digit m gives u = 2.5 + 0.1 m.

if nargin < 1 || isempty(caseName), caseName = 'optimisation'; end
if nargin < 2, overrides = []; end
if nargin < 3 || isempty(studentId), studentId = '52104479'; end

p.case = caseName;
p.student_id = studentId;

% unit factors (display -> SI)
p.cm = 1e-2;  p.mm = 1e-3;  p.MPa = 1e6;  p.GPa = 1e9;

% material and gear constants
p.q = 2.54;  p.kv = 2.1;  p.E = 200*p.GPa;  p.nu = 0.3;
% The brief prints cp = E/(2 pi (1+nu^2)) but only 1-nu^2 reproduces the
% validation sigma_c (see validate_speed_reducer).
p.cp_form = 'one_minus_nu2';

% loading
nDigit = str2double(studentId(end));
mDigit = str2double(studentId(end-1));
switch caseName
    case 'optimisation',   p.T = 1500 + 100*nDigit;  p.u = 2.5 + 0.1*mDigit;
    case 'validation1000', p.T = 1000;  p.u = 3;
    case 'validation2400', p.T = 2400;  p.u = 3;
    case 'reference',      p.T = 1000;  p.u = 3;
    otherwise
        error('speed_reducer_params:badCase', 'Unknown case "%s".', caseName);
end

% mechanical limits (SI)
p.sigma_b_max = 650*p.MPa;  p.sigma_c_max = 800*p.MPa;
p.sigma_s_max = 550*p.MPa;  p.y_max = 0.075*p.mm;

% bounds, display units, order [b m z l1 l2 d1 d2]
p.lb = [2.6; 7.0; 17; 7.3; 7.3; 2.8; 5.0];
p.ub = [4.4; 8.0; 28; 8.3; 8.3; 3.9; 5.5];
p.var_names = {'b','m','z','l1','l2','d1','d2'};
p.var_units = {'cm','mm','-','cm','cm','cm','cm'};
p.z_index = 3;

% gear cross-section constants (Figure 2 of the brief)
p.tooth_offset = 2.4;      % d_st = m (z - 2.4)
p.rim_offset = 4;          % d_w = d_st - 4 m
p.hub_diam_factor = 2.4;   % d_p = 2.4 d
p.web_fraction = 1/3;      % web thickness b/3
p.hub_len_factor = 2;      % hub length 2 d

% geometric requirements
p.l1_min_factor = 1.5;  p.l1_min_offset = 1.9*p.cm;
p.l2_min_factor = 1.1;  p.l2_min_offset = 1.9*p.cm;
p.bm_min = 5;  p.bm_max = 12;
p.overall_max = 160*p.cm;

% solver scaling: objective divided by 1000 cm^3
p.f_scale = 1000*p.cm^3;

if ~isempty(overrides)
    f = fieldnames(overrides);
    for k = 1:numel(f)
        if ~isfield(p, f{k})
            error('speed_reducer_params:badOverride', 'Unknown field "%s".', f{k});
        end
        p.(f{k}) = overrides.(f{k});
    end
end
end
