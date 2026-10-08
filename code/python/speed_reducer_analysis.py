"""speed_reducer_analysis.py

Author: AWD Labs
Student ID: 52104479

Mechanical analysis and volume of the single stage speed reducer.

Unit convention
---------------
Input x (display units): [b cm, m mm, z, l1 cm, l2 cm, d1 cm, d2 cm].
Everything inside is SI: lengths m, stress Pa, torque N m, volume m^3.
The returned dict holds SI values plus a 'display' sub dict
(volume cm^3, stress MPa, deflection mm).

The core function is written without abs()/max() so it also accepts
complex input, which the solver wrapper uses for complex step gradients.
"""
import numpy as np


def _check_inputs(x, p):
    x = np.asarray(x)
    if x.shape != (7,):
        raise ValueError("speed_reducer_analysis:badSize x must have 7 elements, got shape %s" % (x.shape,))
    if not np.all(np.isfinite(x)):
        raise ValueError("speed_reducer_analysis:nonFinite x contains NaN or Inf")
    if not np.all(np.isreal(x)):
        raise ValueError("speed_reducer_analysis:complex x must be real")
    if np.any(x <= 0):
        raise ValueError("speed_reducer_analysis:nonPositive all design variables must be positive")
    for f in ("T", "u", "q", "kv", "E", "nu", "cm", "mm"):
        if f not in p:
            raise KeyError("speed_reducer_analysis:badParams params missing field '%s'" % f)


def speed_reducer_core(x, p):
    """Unchecked analysis, accepts complex x. Returns dict in SI."""
    cm, mm = p["cm"], p["mm"]
    b = x[0] * cm
    m = x[1] * mm
    z = x[2]
    l = np.array([x[3] * cm, x[4] * cm], dtype=x.dtype)
    d = np.array([x[5] * cm, x[6] * cm], dtype=x.dtype)
    T, u = p["T"], p["u"]

    # gear stresses
    sigma_b = p["q"] * 2 * T / (b * m ** 2 * z)
    if p["cp_form"] == "one_minus_nu2":
        cp = p["E"] / (2 * np.pi * (1 - p["nu"] ** 2))
    elif p["cp_form"] == "one_plus_nu2":
        cp = p["E"] / (2 * np.pi * (1 + p["nu"] ** 2))
    else:
        raise ValueError("speed_reducer_analysis:badCp unknown cp_form")
    sigma_c = np.sqrt(cp * p["kv"] * T / (b * m ** 2 * z ** 2) * (1 + u) / u)

    # shafts
    y = 8 / (3 * np.pi) * T * l ** 3 / (p["E"] * d ** 4 * m * z)
    sigma_sb = 16 / np.pi * T * l / (d ** 3 * m * z)
    Tt = np.array([T, u * T])
    sigma_st = 16 / np.pi * Tt / d ** 3
    sigma_s = np.sqrt(sigma_sb ** 2 + 3 * sigma_st ** 2)

    # volumes
    zi = np.array([z, u * z], dtype=x.dtype)
    d_st = m * (zi - p["tooth_offset"])
    d_w = d_st - p["rim_offset"] * m
    d_p = p["hub_diam_factor"] * d
    V_gear = (np.pi * b / 4 * ((d_st ** 2 - d_w ** 2) + p["web_fraction"] * (d_w ** 2 - d_p ** 2))
              + np.pi / 4 * p["hub_len_factor"] * d * (d_p ** 2 - d ** 2))
    V_shaft = np.pi / 4 * d ** 2 * l
    V_total = np.sum(V_gear) + np.sum(V_shaft)

    return dict(V_total=V_total, V_gear=V_gear, V_shaft=V_shaft,
                sigma_b=sigma_b, sigma_c=sigma_c, sigma_s=sigma_s,
                sigma_sb=sigma_sb, sigma_st=sigma_st, y=y,
                d_st=d_st, d_w=d_w, d_p=d_p,
                b=b, m=m, z=z, l=l, d=d,
                b_over_m=b / m, overall=m * z * (1 + u))


def speed_reducer_analysis(x, p):
    """Checked analysis. x in display units, returns dict (SI + display)."""
    _check_inputs(x, p)
    a = speed_reducer_core(np.asarray(x, dtype=float), p)
    a["display"] = dict(V_cm3=a["V_total"] / p["cm"] ** 3,
                        sigma_b_MPa=a["sigma_b"] / p["MPa"],
                        sigma_c_MPa=a["sigma_c"] / p["MPa"],
                        sigma_s_MPa=a["sigma_s"] / p["MPa"],
                        y_mm=a["y"] / p["mm"])
    return a
