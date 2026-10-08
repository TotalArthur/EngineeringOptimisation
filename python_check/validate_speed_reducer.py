"""validate_speed_reducer.py

Author: AWD Labs
Student ID: 52104479

Validation against the brief's data. Runs at T = 1000 and T = 2400 Nm,
and tests the competing readings of the contact stress formula.
Writes results/python/validation.csv.
"""
import os, csv
import numpy as np
from speed_reducer_params import speed_reducer_params
from speed_reducer_analysis import speed_reducer_analysis

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

X_VAL = np.array([3.5, 7.0, 22.0, 7.4, 7.8, 3.5, 5.2])  # b m z l1 l2 d1 d2
TARGETS = [("Volume (cm^3)", 4147.0), ("sigma_b (MPa)", 323.0), ("sigma_c (MPa)", 532.0),
           ("sigma_s1 (MPa)", 512.0), ("sigma_s2 (MPa)", 454.0),
           ("y1 (mm)", 0.0179), ("y2 (mm)", 0.0043)]


def computed(case, overrides=None):
    p = speed_reducer_params(case, overrides)
    a = speed_reducer_analysis(X_VAL, p)
    dsp = a["display"]
    return [dsp["V_cm3"], dsp["sigma_b_MPa"], dsp["sigma_c_MPa"], dsp["sigma_s_MPa"][0],
            dsp["sigma_s_MPa"][1], dsp["y_mm"][0], dsp["y_mm"][1]]


def main():
    os.makedirs(RES, exist_ok=True)
    rows = []
    out = {c: computed(c) for c in ("validation1000", "validation2400")}
    print("%-16s %10s | %10s %8s | %10s %8s" % ("quantity", "target", "T=1000", "err %", "T=2400", "err %"))
    for i, (name, tgt) in enumerate(TARGETS):
        c1, c2 = out["validation1000"][i], out["validation2400"][i]
        e1, e2 = 100 * (c1 - tgt) / tgt, 100 * (c2 - tgt) / tgt
        print("%-16s %10.4g | %10.4g %8.2f | %10.4g %8.2f" % (name, tgt, c1, e1, c2, e2))
        rows.append([name, tgt, c1, e1, c2, e2])
    with open(os.path.join(RES, "validation.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "target", "T1000_computed", "T1000_err_pct", "T2400_computed", "T2400_err_pct"])
        w.writerows(rows)

    print("\nSigma_c readings at T = 2400 Nm (target 532 MPa):")
    cp_rows = []
    for form in ("one_minus_nu2", "one_plus_nu2"):
        sc = computed("validation2400", {"cp_form": form})[2]
        print("  cp = E/(2 pi (%s)): sigma_c = %.2f MPa  (err %.2f %%)" % (form, sc, 100 * (sc - 532) / 532))
        cp_rows.append([form, sc, 100 * (sc - 532) / 532])
    with open(os.path.join(RES, "validation_cp_readings.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cp_form", "sigma_c_MPa", "err_pct_vs_532"])
        w.writerows(cp_rows)
    return out


if __name__ == "__main__":
    main()
