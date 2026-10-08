"""run_optimisation.py

Author: AWD Labs
Student ID: 52104479

Main driver (Python port of run_optimisation.m). Runs, in order:
  1. validation check
  2. multistart continuous solves for both methods (rng seed 1)
  3. integer teeth study z = 17..28
  4. reference case (u = 3, T = 1000 Nm)
Saves .mat and .csv into results/python. Plotting and post processing are
separate scripts (make_figures.py, postprocess.py).
"""
import os, csv, json, sys
import numpy as np
from scipy.io import savemat
from speed_reducer_params import speed_reducer_params
import speed_reducer_problem as prob
import solvers
from validate_speed_reducer import main as validate

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

SEED = 1
N_STARTS = 100        # per method, continuous multistart
N_STARTS_Z = 20       # per method and per z, integer study
N_STARTS_REF = 50


def multistart(p, n_starts, seed, fixed=None, methods=solvers.METHODS):
    rng = np.random.default_rng(seed)
    nfree = len(solvers.free_indices(fixed))
    starts = rng.random((n_starts, nfree))   # uniform in the scaled box
    out = {m: [] for m in methods}
    for k in range(n_starts):
        for m in methods:
            try:
                r = solvers.solve(m, starts[k], p, fixed)
            except Exception as e:  # keep the study going, record the failure
                r = dict(method=m, xs=np.full(nfree, np.nan), x=np.full(7, np.nan), f=np.nan, viol=np.inf,
                         feasible=False, success=False, ok=False, nit=0, nfev=0, time=0.0,
                         history=[], message="exception: %s" % e)
            r["start"] = starts[k]
            out[m].append(r)
    return starts, out


def best_of(runs):
    ok = [r for r in runs if r["ok"]]
    return min(ok, key=lambda r: r["f"]) if ok else None


def write_runs_csv(path, runs):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "ok", "success_flag", "f_cm3", "max_viol", "iterations", "func_evals", "time_s"]
                   + prob_names() + ["message"])
        for i, r in enumerate(runs):
            w.writerow([i, int(r["ok"]), int(r["success"]), r["f"], r["viol"], r["nit"], r["nfev"], r["time"]]
                       + list(r["x"]) + [r["message"]])


def prob_names():
    return ["b_cm", "m_mm", "z", "l1_cm", "l2_cm", "d1_cm", "d2_cm"]


def stats(runs):
    ok = [r for r in runs if r["ok"]]
    f = np.array([r["f"] for r in ok])
    s = dict(n=len(runs), n_ok=len(ok), success_rate=len(ok) / len(runs))
    feas = [r for r in runs if r["feasible"] and np.isfinite(r["f"])]
    if feas:
        fb = min(r["f"] for r in feas)
        s["n_feasible_at_best"] = int(sum(r["f"] <= fb * (1 + 1e-6) for r in feas))
        s["n_feasible"] = len(feas)
    if ok:
        s.update(best=f.min(), worst=f.max(), spread=f.max() - f.min(), mean=f.mean(), std=f.std(ddof=1) if len(f) > 1 else 0.0,
                 n_within_tol=int(np.sum(f <= f.min() * (1 + 1e-4))),
                 mean_iter=np.mean([r["nit"] for r in runs]), mean_nfev=np.mean([r["nfev"] for r in runs]),
                 mean_time=np.mean([r["time"] for r in runs]), total_time=np.sum([r["time"] for r in runs]),
                 median_iter=float(np.median([r["nit"] for r in runs])))
    return s


def main():
    os.makedirs(RES, exist_ok=True)
    print("=== 1. validation ===")
    validate()

    p = speed_reducer_params("optimisation")
    print("\nOptimisation case: T = %g Nm, u = %g (student ID %s)" % (p["T"], p["u"], p["student_id"]))
    summary = {"T": p["T"], "u": p["u"], "options": solvers.OPTIONS, "seed": SEED, "feas_tol": solvers.FEAS_TOL}

    print("\n=== 2. multistart (%d starts per method) ===" % N_STARTS)
    starts, ms = multistart(p, N_STARTS, SEED)
    summary["multistart"] = {}
    best = {}
    for m in solvers.METHODS:
        s = stats(ms[m])
        summary["multistart"][m] = s
        b = best_of(ms[m])
        best[m] = b
        print(m, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in s.items()})
        print("  best x:", np.round(b["x"], 5), "f =", b["f"])
        write_runs_csv(os.path.join(RES, "multistart_%s.csv" % m.replace("-", "_")), ms[m])
    # convergence histories: rerun best start with recording
    hist = {}
    for m in solvers.METHODS:
        runs = ms[m]
        k = next(i for i, r in enumerate(runs) if r is best[m])
        r = solvers.solve(m, runs[k]["start"], p, record=True)
        hist[m] = np.array(r["history"])
        print(m, "history length", len(r["history"]), "re-run f", r["f"])
    # histories from a typical run (median-start) too: first start
    hist_first = {}
    for m in solvers.METHODS:
        r = solvers.solve(m, starts[0], p, record=True)
        hist_first[m] = np.array(r["history"])

    savemat(os.path.join(RES, "multistart.mat"), {
        "starts_scaled": starts,
        **{("f_" + m.replace("-", "_")): np.array([r["f"] for r in ms[m]]) for m in solvers.METHODS},
        **{("ok_" + m.replace("-", "_")): np.array([r["ok"] for r in ms[m]]) for m in solvers.METHODS},
        **{("X_" + m.replace("-", "_")): np.array([r["x"] for r in ms[m]]) for m in solvers.METHODS},
        **{("hist_best_" + m.replace("-", "_")): hist[m] for m in solvers.METHODS},
        **{("hist_start1_" + m.replace("-", "_")): hist_first[m] for m in solvers.METHODS},
        "best_x_sqp": best["sqp"]["x"], "best_x_ip": best["interior-point"]["x"]})

    print("\n=== 3. integer teeth study ===")
    zs = list(range(17, 29))
    zrows = []
    for z in zs:
        _, msz = multistart(p, N_STARTS_Z, SEED, fixed={2: float(z)})
        for m in solvers.METHODS:
            b = best_of(msz[m])
            zrows.append(dict(z=z, method=m, f=(b["f"] if b else np.nan), x=(b["x"] if b else np.full(7, np.nan)),
                              n_ok=sum(r["ok"] for r in msz[m]), n=N_STARTS_Z))
        print("z=%d" % z, {r["method"]: round(r["f"], 4) for r in zrows[-2:]})
    with open(os.path.join(RES, "integer_z.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["z", "method", "V_cm3", "n_feasible_runs", "n_runs"] + prob_names())
        for r in zrows:
            w.writerow([r["z"], r["method"], r["f"], r["n_ok"], r["n"]] + list(r["x"]))
    savemat(os.path.join(RES, "integer_z.mat"), {
        "z": np.array(zs),
        "V_sqp": np.array([r["f"] for r in zrows if r["method"] == "sqp"]),
        "V_ip": np.array([r["f"] for r in zrows if r["method"] == "interior-point"])})

    print("\n=== 4. reference case u = 3, T = 1000 Nm ===")
    pr = speed_reducer_params("reference")
    _, msr = multistart(pr, N_STARTS_REF, SEED)
    summary["reference"] = {}
    for m in solvers.METHODS:
        s = stats(msr[m])
        b = best_of(msr[m])
        summary["reference"][m] = dict(stats=s, best_x=list(b["x"]), best_f=b["f"])
        print(m, "best f =", b["f"], np.round(b["x"], 5), "success", s["success_rate"])
        write_runs_csv(os.path.join(RES, "reference_%s.csv" % m.replace("-", "_")), msr[m])

    summary["best"] = {m: dict(x=list(best[m]["x"]), f=best[m]["f"]) for m in solvers.METHODS}
    with open(os.path.join(RES, "run_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print("\nsaved to", RES)


if __name__ == "__main__":
    main()
