function a = speed_reducer_core(x, p)
% SPEED_REDUCER_CORE  Unchecked analysis in SI. Accepts complex x (complex step).
%
% Author: AWD Labs
% Student ID: 52104479
%
% Unit convention: x in display units [b cm, m mm, z, l1 cm, l2 cm, d1 cm, d2 cm];
% everything computed in SI. No abs()/max() so complex input is safe.

cm = p.cm;  mm = p.mm;
b = x(1)*cm;  m = x(2)*mm;  z = x(3);
l = [x(4); x(5)]*cm;
d = [x(6); x(7)]*cm;
T = p.T;  u = p.u;

% gear stresses
sigma_b = p.q*2*T/(b*m^2*z);
switch p.cp_form
    case 'one_minus_nu2', cp = p.E/(2*pi*(1 - p.nu^2));
    case 'one_plus_nu2',  cp = p.E/(2*pi*(1 + p.nu^2));
    otherwise, error('speed_reducer_core:badCp', 'Unknown cp_form.');
end
sigma_c = sqrt(cp*p.kv*T/(b*m^2*z^2)*(1 + u)/u);

% shafts
y = 8/(3*pi)*T*l.^3./(p.E*d.^4*m*z);
sigma_sb = 16/pi*T*l./(d.^3*m*z);
sigma_st = 16/pi*[T; u*T]./d.^3;
sigma_s = sqrt(sigma_sb.^2 + 3*sigma_st.^2);

% volumes
zi = [z; u*z];
d_st = m*(zi - p.tooth_offset);
d_w = d_st - p.rim_offset*m;
d_p = p.hub_diam_factor*d;
V_gear = pi*b/4*((d_st.^2 - d_w.^2) + p.web_fraction*(d_w.^2 - d_p.^2)) ...
       + pi/4*p.hub_len_factor*d.*(d_p.^2 - d.^2);
V_shaft = pi/4*d.^2.*l;

a.V_total = sum(V_gear) + sum(V_shaft);
a.V_gear = V_gear;  a.V_shaft = V_shaft;
a.sigma_b = sigma_b;  a.sigma_c = sigma_c;  a.sigma_s = sigma_s;
a.sigma_sb = sigma_sb;  a.sigma_st = sigma_st;  a.y = y;
a.d_st = d_st;  a.d_w = d_w;  a.d_p = d_p;
a.b = b;  a.m = m;  a.z = z;  a.l = l;  a.d = d;
a.b_over_m = b/m;
a.overall = m*z*(1 + u);
end
