"""speed_reducer_problem.py

Author: AWD Labs
Student ID: 52104479

Objective, normalised constraints, scaling and complex step derivatives.
Constraints are g_i = value/limit - 1 <= 0 so every row is O(1). The linear
geometric rows use a constant reference scale instead (see constraint_values).

Scaled variables: xs = (x - lb) / (ub - lb), so the box becomes [0, 1]^n.
"""
import numpy as np
from speed_reducer_analysis import speed_reducer_core

CON_NAMES = [
    "gear bending stress", "gear contact stress",
    "shaft 1 combined stress", "shaft 2 combined stress",
    "shaft 1 deflection", "shaft 2 deflection",
    "shaft 1 min length", "shaft 2 min length",
    "b/m lower limit", "b/m upper limit", "overall size",
]


def to_scaled(x, p):
    return (np.asarray(x) - p["lb"]) / (p["ub"] - p["lb"])


def from_scaled(xs, p):
    return p["lb"] + xs * (p["ub"] - p["lb"])


def _fixed_full(xs_red, p, fixed):
    """Insert fixed variables. fixed = dict {index: display value} or None."""
    if not fixed:
        return from_scaled(xs_red, p)
    n = 7
    free = [i for i in range(n) if i not in fixed]
    xs = np.zeros(n, dtype=np.asarray(xs_red).dtype)
    xs[free] = xs_red
    x = from_scaled(xs, p)
    for i, v in fixed.items():
        x = x.astype(xs.dtype) if x.dtype != xs.dtype else x
        x[i] = v
    return x


def objective_value(xs, p, fixed=None):
    x = _fixed_full(xs, p, fixed)
    a = speed_reducer_core(x, p)
    return a["V_total"] / p["f_scale"]


def constraint_values(xs, p, fixed=None):
    """Normalised inequality constraints g <= 0 (11 rows)."""
    x = _fixed_full(xs, p, fixed)
    a = speed_reducer_core(x, p)
    cm = p["cm"]
    d, l = a["d"], a["l"]
    # The four geometric rows are linear in x, so they are normalised by a
    # constant reference (not by a variable) to keep them linear and convex.
    ref_l = p["lb"][3] * cm          # smallest allowed shaft length
    ref_b = p["lb"][0] * cm          # smallest allowed gear width
    m_, b_ = a["m"], a["b"]
    g = [a["sigma_b"] / p["sigma_b_max"] - 1,
         a["sigma_c"] / p["sigma_c_max"] - 1,
         a["sigma_s"][0] / p["sigma_s_max"] - 1,
         a["sigma_s"][1] / p["sigma_s_max"] - 1,
         a["y"][0] / p["y_max"] - 1,
         a["y"][1] / p["y_max"] - 1,
         (p["l1_min_factor"] * d[0] + p["l1_min_offset"] - l[0]) / ref_l,
         (p["l2_min_factor"] * d[1] + p["l2_min_offset"] - l[1]) / ref_l,
         (p["bm_min"] * m_ - b_) / ref_b,
         (b_ - p["bm_max"] * m_) / ref_b,
         a["overall"] / p["overall_max"] - 1]
    return np.array(g, dtype=np.asarray(x).dtype)


def cs_gradient(fun, xs, h=1e-30):
    """Complex step Jacobian of fun (scalar or vector valued) at real xs."""
    xs = np.asarray(xs, dtype=float)
    f0 = np.atleast_1d(fun(xs.astype(complex)))
    J = np.zeros((f0.size, xs.size))
    for k in range(xs.size):
        xc = xs.astype(complex)
        xc[k] += 1j * h
        J[:, k] = np.imag(np.atleast_1d(fun(xc))) / h
    return J


def objective_grad(xs, p, fixed=None):
    return cs_gradient(lambda v: objective_value(v, p, fixed), xs)[0]


def constraint_jac(xs, p, fixed=None):
    return cs_gradient(lambda v: constraint_values(v, p, fixed), xs)


def max_violation(xs, p, fixed=None):
    g = constraint_values(np.asarray(xs, float), p, fixed)
    return max(0.0, float(np.max(g)))
