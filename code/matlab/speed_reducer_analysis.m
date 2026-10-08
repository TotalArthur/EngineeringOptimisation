function a = speed_reducer_analysis(x, p)
% SPEED_REDUCER_ANALYSIS  Checked analysis: volume, stresses, deflections, geometry.
%
% Author: AWD Labs
% Student ID: 52104479
%
% UNIT CONVENTION
%   Input  x (display units): [b cm, m mm, z, l1 cm, l2 cm, d1 cm, d2 cm]
%   Inside: SI (m, Pa, N m, m^3).
%   Output: struct a in SI plus a.display (volume cm^3, stress MPa, deflection mm).
%
% Errors: speed_reducer_analysis:badSize, :nonFinite, :nonPositive, :badParams

x = x(:);
if numel(x) ~= 7
    error('speed_reducer_analysis:badSize', 'x must have 7 elements, got %d.', numel(x));
end
if ~isnumeric(x) || ~isreal(x) || any(~isfinite(x))
    error('speed_reducer_analysis:nonFinite', 'x must be real and finite.');
end
if any(x <= 0)
    error('speed_reducer_analysis:nonPositive', 'All design variables must be positive.');
end
need = {'T','u','q','kv','E','nu','cm','mm','cp_form'};
for k = 1:numel(need)
    if ~isfield(p, need{k})
        error('speed_reducer_analysis:badParams', 'params missing field "%s".', need{k});
    end
end

a = speed_reducer_core(x, p);
a.display.V_cm3 = a.V_total/p.cm^3;
a.display.sigma_b_MPa = a.sigma_b/p.MPa;
a.display.sigma_c_MPa = a.sigma_c/p.MPa;
a.display.sigma_s_MPa = a.sigma_s/p.MPa;
a.display.y_mm = a.y/p.mm;
end
