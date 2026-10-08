"""postprocess.py

Author: AWD Labs
Student ID: 52104479

Post processing of the best design per method: constraint and bound
tables, Lagrange multipliers, first order optimality, Hessian of the
Lagrangian, sampled convexity test and an independent sanity check of the
original (unnormalised) constraints. Reads results/python/run_summary.json.
Writes results/python/postprocess.json and several csv files.
"""
import os, json, csv
import numpy as np
from scipy.optimize import nnls
from scipy.io import savemat
from speed_reducer_params import speed_reducer_params
from speed_reducer_analysis import speed_reducer_analysis
import speed_reducer_problem as prob
import solvers

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "results", "python")
ACTIVE_TOL = 1e-4      # relative physical margin below this counts as active
BOUND_TOL = 1e-4       # distance (fraction of range) counted as "at bound"
# (interior-point stops about 1e-5 inside the feasible set, so 1e-4 is used for both methods)
N_SEG = 5000
CON_LIMITS = None


def physical_constraints(x, p):
    """Original, unnormalised constraints written out directly from the analysis
    (independent of the normalised g). Returns list of (name, value, limit, unit, sense)."""
    a = speed_reducer_analysis(x, p)
    d = a["display"]
    cm = p["cm"]
    return [
        ("gear bending stress", d["sigma_b_MPa"], p["sigma_b_max"] / p["MPa"], "MPa", "<="),
        ("gear contact stress", d["sigma_c_MPa"], p["sigma_c_max"] / p["MPa"], "MPa", "<="),
        ("shaft 1 combined stress", d["sigma_s_MPa"][0], p["sigma_s_max"] / p["MPa"], "MPa", "<="),
        ("shaft 2 combined stress", d["sigma_s_MPa"][1], p["sigma_s_max"] / p["MPa"], "MPa", "<="),
        ("shaft 1 deflection", d["y_mm"][0], p["y_max"] / p["mm"], "mm", "<="),
        ("shaft 2 deflection", d["y_mm"][1], p["y_max"] / p["mm"], "mm", "<="),
        ("shaft 1 min length (1.5 d1 + 1.9)", p["l1_min_factor"] * x[5] + p["l1_min_offset"] / cm, x[3], "cm", "<="),
        ("shaft 2 min length (1.1 d2 + 1.9)", p["l2_min_factor"] * x[6] + p["l2_min_offset"] / cm, x[4], "cm", "<="),
        ("b/m lower limit", p["bm_min"], a["b_over_m"], "-", "<="),
        ("b/m upper limit", a["b_over_m"], p["bm_max"], "-", "<="),
        ("overall size m z (1+u)", a["overall"] / cm, p["overall_max"] / cm, "cm", "<="),
    ]


def kkt_multipliers(xs, p):
    """Non negative least squares estimate of the multipliers at xs (scaled space).
    Active set: constraints with relative margin < ACTIVE_TOL, bounds within BOUND_TOL."""
    g = prob.constraint_values(xs, p)
    Jg = prob.constraint_jac(xs, p)
    gf = prob.objective_grad(xs, p)
    x = prob.from_scaled(xs, p)
    rel = [(lim - val) / abs(lim) for (_, val, lim, _, _) in physical_constraints(x, p)]
    act_c = [i for i in range(len(g)) if rel[i] < ACTIVE_TOL]
    act_lo = [i for i in range(7) if xs[i] < BOUND_TOL]
    act_up = [i for i in range(7) if xs[i] > 1 - BOUND_TOL]
    cols = [Jg[i] for i in act_c]
    for i in act_lo:
        e = np.zeros(7); e[i] = -1.0; cols.append(e)    # lower bound: -(x_i) <= 0
    for i in act_up:
        e = np.zeros(7); e[i] = 1.0; cols.append(e)     # upper bound: x_i - 1 <= 0
    A = np.array(cols).T                                # 7 x n_active
    lam, _ = nnls(A, -gf)
    resid = gf + A @ lam
    return dict(g=g, Jg=Jg, gf=gf, act_c=act_c, act_lo=act_lo, act_up=act_up, A=A, lam=lam,
                lam_c=lam[:len(act_c)], lam_lo=lam[len(act_c):len(act_c) + len(act_lo)],
                lam_up=lam[len(act_c) + len(act_lo):], first_order=float(np.max(np.abs(resid))),
                resid_norm=float(np.linalg.norm(resid)))


def lagrangian_grad(xs, p, lam_full):
    return prob.objective_grad(xs, p) + prob.constraint_jac(xs, p).T @ lam_full


def hessian_lagrangian(xs, p, lam_full, h=1e-6):
    n = 7
    H = np.zeros((n, n))
    for j in range(n):
        xp, xm = xs.copy(), xs.copy()
        xp[j] += h; xm[j] -= h
        H[:, j] = (lagrangian_grad(xp, p, lam_full) - lagrangian_grad(xm, p, lam_full)) / (2 * h)
    return 0.5 * (H + H.T)


def convexity_test(p, seed=2):
    """Midpoint convexity test along random segments inside the bounds.
    For each function h (objective and 11 constraints): count segments where
    h(mid) > (h(a)+h(b))/2 + tol, i.e. convexity is violated."""
    rng = np.random.default_rng(seed)
    names = ["objective"] + prob.CON_NAMES
    viol = np.zeros(12, dtype=int)
    worst = np.zeros(12)
    for _ in range(N_SEG):
        a, b = rng.random(7), rng.random(7)
        m = 0.5 * (a + b)
        fa = np.r_[prob.objective_value(a, p), prob.constraint_values(a, p)]
        fb = np.r_[prob.objective_value(b, p), prob.constraint_values(b, p)]
        fm = np.r_[prob.objective_value(m, p), prob.constraint_values(m, p)]
        gap = fm - 0.5 * (fa + fb)       # > 0 means non convex along this segment
        scale = np.maximum(1e-12, 0.5 * (np.abs(fa) + np.abs(fb)))
        bad = gap > 1e-9 * scale
        viol += bad
        worst = np.maximum(worst, gap / scale)
    return names, viol, worst


def analyse(label, p, x, method, rows_out):
    xs = prob.to_scaled(x, p)
    k = kkt_multipliers(xs, p)
    lam_full = np.zeros(11)
    for idx, l in zip(k["act_c"], k["lam_c"]):
        lam_full[idx] = l
    H = hessian_lagrangian(xs, p, lam_full)
    eigH = np.linalg.eigvalsh(H)
    # reduced Hessian on null space of active constraint + bound gradients
    A = k["A"]
    if A.shape[1] > 0:
        U, S, Vt = np.linalg.svd(A.T)
        rank = int(np.sum(S > 1e-9))
        Z = Vt[rank:].T
    else:
        Z = np.eye(7)
    Hred = Z.T @ H @ Z if Z.shape[1] else np.zeros((0, 0))
    eigZ = np.linalg.eigvalsh(Hred) if Z.shape[1] else np.array([])

    a = speed_reducer_analysis(x, p)
    out = dict(label=label, method=method, x=list(x), V_cm3=a["display"]["V_cm3"])
    # constraint table
    g = k["g"]
    phys = physical_constraints(x, p)
    ctab = []
    for i, (nm, val, lim, unit, sense) in enumerate(phys):
        margin = 100 * (lim - val) / abs(lim)
        ctab.append(dict(name=nm, value=val, limit=lim, unit=unit, g=float(g[i]), margin_pct=float(margin),
                         active=bool(i in k["act_c"]), satisfied=bool(val <= lim * (1 + 1e-9) + 1e-12),
                         multiplier=float(lam_full[i])))
    out["constraints"] = ctab
    out["all_original_constraints_satisfied"] = all(c["satisfied"] for c in ctab)
    # bounds table
    btab = []
    for i in range(7):
        status = "lower" if xs[i] < BOUND_TOL else ("upper" if xs[i] > 1 - BOUND_TOL else "free")
        mult = 0.0
        if i in k["act_lo"]:
            mult = float(k["lam_lo"][k["act_lo"].index(i)])
        if i in k["act_up"]:
            mult = float(k["lam_up"][k["act_up"].index(i)])
        btab.append(dict(name=p["var_names"][i], unit=p["var_units"][i], value=float(x[i]), lb=float(p["lb"][i]),
                         ub=float(p["ub"][i]), status=status, multiplier_scaled=mult))
    out["bounds"] = btab
    out["first_order_optimality"] = k["first_order"]
    out["n_active_constraints"] = len(k["act_c"])
    out["n_active_bounds"] = len(k["act_lo"]) + len(k["act_up"])
    out["eig_hessian_lagrangian"] = list(map(float, eigH))
    out["eig_reduced_hessian"] = list(map(float, eigZ))
    out["null_space_dim"] = int(Z.shape[1])
    out["lagrangian_scaling_note"] = "objective scaled by 1000 cm^3, variables scaled to [0,1]"
    rows_out.append(out)
    return out


def main():
    summary = json.load(open(os.path.join(RES, "run_summary.json")))
    p = speed_reducer_params("optimisation")
    pr = speed_reducer_params("reference")
    results = []
    for m in solvers.METHODS:
        analyse("optimisation", p, np.array(summary["best"][m]["x"]), m, results)
    for m in solvers.METHODS:
        analyse("reference", pr, np.array(summary["reference"][m]["best_x"]), m, results)

    names, viol, worst = convexity_test(p)
    conv = [dict(function=n, n_violations=int(v), n_segments=N_SEG, frac=float(v) / N_SEG, worst_rel_gap=float(w))
            for n, v, w in zip(names, viol, worst)]
    convref = None

    # sensitivity to the printed (1+nu^2) cp reading: does the optimum move?
    p2 = speed_reducer_params("optimisation", {"cp_form": "one_plus_nu2"})
    r2 = solvers.solve("sqp", np.full(7, 0.5), p2)
    sens = dict(cp_form="one_plus_nu2", f=r2["f"], x=list(r2["x"]), feasible=bool(r2["feasible"]),
                sigma_c_MPa=float(speed_reducer_analysis(r2["x"], p2)["display"]["sigma_c_MPa"]))

    with open(os.path.join(RES, "postprocess.json"), "w") as f:
        json.dump(dict(results=results, convexity=conv, cp_sensitivity=sens, active_tol=ACTIVE_TOL,
                       bound_tol=BOUND_TOL), f, indent=2)
    for tag in ("optimisation", "reference"):
        for m in solvers.METHODS:
            r = next(x for x in results if x["label"] == tag and x["method"] == m)
            base = "%s_%s" % (tag, m.replace("-", "_"))
            with open(os.path.join(RES, "constraints_%s.csv" % base), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["constraint", "value", "limit", "unit", "g_normalised", "margin_pct", "active", "multiplier"])
                for c in r["constraints"]:
                    w.writerow([c["name"], c["value"], c["limit"], c["unit"], c["g"], c["margin_pct"], int(c["active"]), c["multiplier"]])
            with open(os.path.join(RES, "bounds_%s.csv" % base), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["variable", "unit", "value", "lb", "ub", "status", "multiplier_scaled"])
                for b in r["bounds"]:
                    w.writerow([b["name"], b["unit"], b["value"], b["lb"], b["ub"], b["status"], b["multiplier_scaled"]])
    with open(os.path.join(RES, "convexity_samples.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function", "n_violations", "n_segments", "fraction", "worst_relative_gap"])
        for c in conv:
            w.writerow([c["function"], c["n_violations"], c["n_segments"], c["frac"], c["worst_rel_gap"]])
    savemat(os.path.join(RES, "postprocess.mat"), {
        "eig_H_sqp": results[0]["eig_hessian_lagrangian"], "eig_H_ip": results[1]["eig_hessian_lagrangian"],
        "eig_Hred_sqp": results[0]["eig_reduced_hessian"], "eig_Hred_ip": results[1]["eig_reduced_hessian"]})

    # print
    for r in results:
        print("\n==== %s / %s : V = %.4f cm^3, first-order optimality %.2e ====" % (r["label"], r["method"], r["V_cm3"], r["first_order_optimality"]))
        for c in r["constraints"]:
            print("  %-36s %10.5g %s %10.5g  margin %8.3f %%  %s  lam=%.4g" % (c["name"], c["value"], c["unit"], c["limit"], c["margin_pct"], "ACTIVE" if c["active"] else "inactive", c["multiplier"]))
        for b in r["bounds"]:
            print("  %-3s %9.5f %-3s [%g, %g] %s lam=%.4g" % (b["name"], b["value"], b["unit"], b["lb"], b["ub"], b["status"], b["multiplier_scaled"]))
        print("  eig H_L:", np.round(r["eig_hessian_lagrangian"], 5))
        print("  eig reduced H (dim %d):" % r["null_space_dim"], np.round(r["eig_reduced_hessian"], 5))
        print("  all original constraints satisfied:", r["all_original_constraints_satisfied"])
    print("\nConvexity (midpoint test, %d segments)" % N_SEG)
    for c in conv:
        print("  %-28s violations %5d (%.1f %%)  worst rel gap %.3g" % (c["function"], c["n_violations"], 100 * c["frac"], c["worst_rel_gap"]))
    print("\nsensitivity to cp reading:", sens)


if __name__ == "__main__":
    main()
