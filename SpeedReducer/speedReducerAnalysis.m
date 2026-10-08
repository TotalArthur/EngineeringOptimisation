function r = speedReducerAnalysis(x, params)
% speedReducerAnalysis  Mechanical analysis and material volume of the speed reducer.
% Author: AWD Labs
% EG503X Engineering Optimisation, Assignment 1
%
% INPUTS
%   x      design vector in the units of the brief:
%          x = [b (cm), m (mm), z (-), l1 (cm), l2 (cm), d1 (cm), d2 (cm)]
%   params struct made at the top of main.m (everything is in SI units)
%
% OUTPUT
%   r      struct with stresses (Pa), deflections (m), volumes (m^3)
%          and the 11 constraint values and limits (SI units)
%
% The solver works with x in cm and mm because then every variable is
% between about 2 and 30, which is much better scaled than metres
% (0.03 next to 22). All the maths below is done in SI.

%% unpack the design variables and convert to SI
xSI = x .* params.xToSI;
b = xSI(1);              % gear width (m)
m = xSI(2);              % gear module (m)
z = xSI(3);              % teeth on gear 1 (-)
l = [xSI(4), xSI(5)];    % shaft lengths between bearings (m)
d = [xSI(6), xSI(7)];    % shaft diameters (m)

%% unpack fixed parameters
T  = params.T;           % input torque on shaft 1 (Nm)
u  = params.u;           % speed ratio
E  = params.E;           % Young's modulus (Pa)
nu = params.nu;          % Poisson's ratio

%% gear tooth stresses
% Bending stress at the tooth root. Tangential force is 2T/(m z), spread
% over the face width b and the module m, times the form factor q.
r.sigma_b = params.q * 2*T / (b * m^2 * z);

% Contact (compressive) stress. cp is the elastic coefficient for two gears
% made of the same material. The (1 - nu^2) version is the standard Hertz
% contact result (Shigley). The sign is a switch in params so the other
% version in the brief can be tested.
cp = E / (2*pi * (1 + params.cpNuSign * nu^2));
r.sigma_c = sqrt( cp * params.kv * T / (b * m^2 * z^2) * (1 + u) / u );

%% shaft deflection and stresses (i = 1, 2 stored as columns 1 and 2)
r.y = 8/(3*pi) * T * l.^3 ./ (E * d.^4 * m * z);

sigma_sb = 16/pi * T * l ./ (d.^3 * m * z);          % bending stress
sigma_st = 16/pi * T * [1, u] ./ d.^3;               % torsional stress (shaft 2 carries u*T)
r.sigma_s = sqrt(sigma_sb.^2 + 3 * sigma_st.^2);     % combined stress

%% volume
% Gear 1 has z teeth, gear 2 has u*z teeth.
zGear = [z, u*z];
d_st = m * (zGear - params.kTip);     % outer (tip) diameter of each gear
d_w  = d_st - params.kRim * m;        % inner diameter of the rim
d_p  = params.kHub * d;               % outer diameter of the hub

% rim + web + hub, one value for each gear
rim_web = pi*b/4 * ( (d_st.^2 - d_w.^2) + params.kWeb * (d_w.^2 - d_p.^2) );
hub     = pi/4 * (params.kHubWidth * d) .* (d_p.^2 - d.^2);
r.V_gears = rim_web + hub;

% shafts treated as plain cylinders
r.V_shafts = pi/4 * d.^2 .* l;

r.V_gear_total  = sum(r.V_gears);
r.V_shaft_total = sum(r.V_shafts);
r.V_total       = r.V_gear_total + r.V_shaft_total;

%% constraints
% Every constraint is stored as a value and a limit. The order here must
% match params.conNames, params.conIsMin and params.conToDisp in main.m.
r.conValue = [r.sigma_b, r.sigma_c, r.sigma_s(1), r.sigma_s(2), ...
              r.y(1), r.y(2), l(1), l(2), b/m, b/m, m*z*(1+u)];

r.conLimit = [params.sigma_b_max, params.sigma_c_max, ...
              params.sigma_s_max, params.sigma_s_max, ...
              params.y_max, params.y_max, ...
              params.l1_perD1 * d(1) + params.l1_extra, ...
              params.l2_perD2 * d(2) + params.l2_extra, ...
              params.bm_min, params.bm_max, params.size_max];
end
