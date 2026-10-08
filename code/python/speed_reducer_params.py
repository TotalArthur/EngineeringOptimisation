"""speed_reducer_params.py

Author: AWD Labs
Student ID: 52104479
Course: EG503X Engineering Optimisation, Assignment 1 (speed reducer)

Python port of speed_reducer_params.m. Returns a dict holding every
constant, limit and bound. No magic numbers live anywhere else.

Unit convention (whole code base)
---------------------------------
Design vector x is given in DISPLAY units:
    x = [b (cm), m (mm), z (-), l1 (cm), l2 (cm), d1 (cm), d2 (cm)]
speed_reducer_analysis converts to SI (m, Pa, N m, m^3) straight away and
works in SI internally. Results are converted back only for display.
"""
import numpy as np

STUDENT_ID = "52104479"

VAR_NAMES = ["b", "m", "z", "l1", "l2", "d1", "d2"]
VAR_UNITS = ["cm", "mm", "-", "cm", "cm", "cm", "cm"]


def speed_reducer_params(case="optimisation", overrides=None, student_id=STUDENT_ID):
    """Build the parameter dict.

    case:
      'optimisation'   T = 1500 + 100 n, u = 2.5 + 0.1 m from the student ID
      'validation1000' brief's validation data, T = 1000 Nm, u = 3
      'validation2400' same data but T = 2400 Nm
      'reference'      u = 3, T = 1000 Nm with the course limits and bounds
    overrides: dict of fields to replace (e.g. {'cp_form': 'one_plus_nu2'}).
    """
    p = {}
    p["case"] = case
    p["student_id"] = student_id

    # unit conversion factors (display -> SI)
    p["cm"] = 1e-2
    p["mm"] = 1e-3
    p["MPa"] = 1e6
    p["GPa"] = 1e9

    # material and gear constants (same for all cases)
    p["q"] = 2.54
    p["kv"] = 2.1
    p["E"] = 200 * p["GPa"]
    p["nu"] = 0.3
    # elastic coefficient form. The brief prints 1+nu^2 but only 1-nu^2
    # reproduces the validation sigma_c (see validate_speed_reducer).
    p["cp_form"] = "one_minus_nu2"

    # loading
    n = int(student_id[-1])
    mm_digit = int(student_id[-2])
    if case == "optimisation":
        p["T"] = 1500.0 + 100.0 * n
        p["u"] = 2.5 + 0.1 * mm_digit
    elif case == "validation1000":
        p["T"] = 1000.0
        p["u"] = 3.0
    elif case == "validation2400":
        p["T"] = 2400.0
        p["u"] = 3.0
    elif case == "reference":
        p["T"] = 1000.0
        p["u"] = 3.0
    else:
        raise ValueError("speed_reducer_params:badCase unknown case '%s'" % case)

    # limits (SI)
    p["sigma_b_max"] = 650 * p["MPa"]
    p["sigma_c_max"] = 800 * p["MPa"]
    p["sigma_s_max"] = 550 * p["MPa"]
    p["y_max"] = 0.075 * p["mm"]

    # bounds in display units, order [b m z l1 l2 d1 d2]
    p["lb"] = np.array([2.6, 7.0, 17.0, 7.3, 7.3, 2.8, 5.0])
    p["ub"] = np.array([4.4, 8.0, 28.0, 8.3, 8.3, 3.9, 5.5])
    p["var_names"] = VAR_NAMES
    p["var_units"] = VAR_UNITS
    p["z_index"] = 2

    # geometry constants (gear cross-section, Figure 2 of the brief)
    p["tooth_offset"] = 2.4      # d_st = m (z - 2.4)
    p["rim_offset"] = 4.0        # d_w = d_st - 4 m
    p["hub_diam_factor"] = 2.4   # d_p = 2.4 d
    p["web_fraction"] = 1.0 / 3.0  # web thickness b/3
    p["hub_len_factor"] = 2.0    # hub length 2 d

    # geometric requirements (lengths in cm as in the brief)
    p["l1_min_factor"] = 1.5
    p["l1_min_offset"] = 1.9 * p["cm"]
    p["l2_min_factor"] = 1.1
    p["l2_min_offset"] = 1.9 * p["cm"]
    p["bm_min"] = 5.0
    p["bm_max"] = 12.0
    p["overall_max"] = 160.0 * p["cm"]

    # solver scaling
    p["f_scale"] = 1000.0 * (p["cm"] ** 3)  # objective divided by 1000 cm^3

    if overrides:
        for k, v in overrides.items():
            if k not in p:
                raise KeyError("speed_reducer_params:badOverride unknown field '%s'" % k)
            p[k] = v
    return p
