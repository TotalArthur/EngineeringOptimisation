"""literature_comparison.py

Author: AWD Labs
Student ID: 52104479

Compares the optima with Golinski (1970) Table 1. Golinski lists module in
cm; converted to mm here. Evaluates his two designs with this code's model
(volume formula of the brief) and reports which of OUR bounds and limits
they break. Writes results/python/literature.csv.
"""
import os, json, csv
import numpy as np
from speed_reducer_params import speed_reducer_params
from speed_reducer_analysis import speed_reducer_analysis
import speed_reducer_problem as prob

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

# [b cm, m mm, z, l1, l2, d1, d2], f reported by Golinski (cm^3)
GOLINSKI = {
    "Golinski crude Monte Carlo": (np.array([4.4, 6.0, 17, 7.3, 8.1, 3.4, 5.0]), 2236.35),
    "Golinski stray process": (np.array([3.6, 7.0, 18, 6.6, 8.2, 2.8, 5.2]), 2247.79),
}


def main():
    summ = json.load(open(os.path.join(RES, "run_summary.json")))
    pr = speed_reducer_params("reference")
    po = speed_reducer_params("optimisation")
    rows = []
    designs = {
        "This work, T=2400 Nm, u=3.2": (np.array(summ["best"]["sqp"]["x"]), None, po),
        "Reference, T=1000 Nm, u=3": (np.array(summ["reference"]["sqp"]["best_x"]), None, pr),
    }
    for k, (x, f) in GOLINSKI.items():
        designs[k] = (x, f, pr)
    for name, (x, fpaper, p) in designs.items():
        a = speed_reducer_analysis(x, p)
        xs = prob.to_scaled(x, p)
        g = prob.constraint_values(xs, p)
        out_of_bounds = [p["var_names"][i] for i in range(7) if x[i] < p["lb"][i] - 1e-9 or x[i] > p["ub"][i] + 1e-9]
        at_bounds = [p["var_names"][i] for i in range(7) if abs(x[i] - p["lb"][i]) < 1e-6 or abs(x[i] - p["ub"][i]) < 1e-6]
        viol = [prob.CON_NAMES[i] for i in range(11) if g[i] > 1e-6]
        active = [prob.CON_NAMES[i] for i in range(11) if abs(g[i]) <= 1e-4]
        rows.append(dict(design=name, f_paper=fpaper, f_model=a["display"]["V_cm3"], x=list(map(float, x)),
                         outside_bounds=out_of_bounds, at_bound=at_bounds, violated=viol, active=active))
        print(name, "f_model=%.2f" % a["display"]["V_cm3"], "f_paper", fpaper)
        print("   x =", np.round(x, 4), "\n   outside our bounds:", out_of_bounds, "\n   at bounds:", at_bounds,
              "\n   violates constraints (course limits, T=1000, u=3):", viol, "\n   active:", active)
    with open(os.path.join(RES, "literature.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "f_reported_cm3", "f_this_model_cm3", "b_cm", "m_mm", "z", "l1_cm", "l2_cm", "d1_cm", "d2_cm",
                    "outside_our_bounds", "at_bounds", "violated_constraints", "active_constraints"])
        for r in rows:
            w.writerow([r["design"], r["f_paper"], r["f_model"]] + r["x"] + [";".join(r["outside_bounds"]),
                       ";".join(r["at_bound"]), ";".join(r["violated"]), ";".join(r["active"])])
    json.dump(rows, open(os.path.join(RES, "literature.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
