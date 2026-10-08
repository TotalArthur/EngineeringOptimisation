"""make_figures.py

Author: AWD Labs
Student ID: 52104479

Report figures (PNG 300 dpi + PDF) from results/python. Kept separate from
the solver. Style: one blue/orange pair (colourblind safe) plus markers so
the two methods never rely on colour alone.
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

C = {"sqp": "#1f5fa8", "interior-point": "#d9731a"}
LBL = {"sqp": "SQP (SLSQP)", "interior-point": "Interior point (trust-constr)"}
MK = {"sqp": "o", "interior-point": "s"}
plt.rcParams.update({"font.size": 11, "axes.labelsize": 11, "axes.titlesize": 11, "legend.fontsize": 9,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
                     "figure.dpi": 100, "pdf.fonttype": 42})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=300)
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    plt.close(fig)


def fig_convergence(ms, post):
    fstar = post["best_sqp"]
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.4))
    for m, key in (("sqp", "sqp"), ("interior-point", "interior_point")):
        h = np.atleast_2d(ms["hist_start1_" + key])
        it = np.arange(1, len(h) + 1)
        ax[0].semilogy(it, np.maximum(np.abs(h[:, 0] - fstar), 1e-9), marker=MK[m], ms=4, lw=1.4, color=C[m], label=LBL[m])
        ax[1].semilogy(it, np.maximum(h[:, 1], 1e-16), marker=MK[m], ms=4, lw=1.4, color=C[m], label=LBL[m])
    ax[0].set_xlabel("Iteration"); ax[0].set_ylabel("|Volume - best known| (cm$^3$)")
    ax[0].set_title("(a) Objective error")
    ax[1].set_xlabel("Iteration"); ax[1].set_ylabel("Max normalised violation (-)")
    ax[1].set_title("(b) Constraint violation")
    ax[0].legend(frameon=False)
    save(fig, "fig1_convergence")


def fig_multistart(ms):
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.4))
    fbest = float(np.min(ms["f_sqp"]))
    for m, key in (("sqp", "sqp"), ("interior-point", "interior_point")):
        f = np.asarray(ms["f_" + key]).ravel()
        ok = np.asarray(ms["ok_" + key]).ravel().astype(bool)
        idx = np.arange(1, len(f) + 1)
        ax[0].plot(idx[ok], f[ok] - fbest, MK[m], ms=3.5, color=C[m], label=LBL[m], alpha=0.85)
        ax[0].plot(idx[~ok], f[~ok] - fbest, "x", ms=6, color="k")
    ax[0].set_xlabel("Start number"); ax[0].set_ylabel("Final volume above best (cm$^3$)")
    ax[0].set_title("(a) Final objective per start"); ax[0].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.92))
    best = ms["f_sqp"].min()
    for m, key in (("sqp", "sqp"), ("interior-point", "interior_point")):
        f = np.asarray(ms["f_" + key]).ravel()
        ax[1].hist(f - best, bins=20, color=C[m], alpha=0.6, label=LBL[m], edgecolor="white")
    ax[1].set_xlabel("Final volume above best (cm$^3$)"); ax[1].set_ylabel("Number of starts")
    ax[1].set_title("(b) Spread of final values")
    save(fig, "fig2_multistart")


def fig_integer_z(iz):
    z = np.asarray(iz["z"]).ravel()
    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    ax.plot(z, iz["V_sqp"].ravel(), "o-", color=C["sqp"], ms=5, lw=1.4, label=LBL["sqp"])
    ax.plot(z, iz["V_ip"].ravel(), "s--", color=C["interior-point"], ms=4, lw=1.2, label=LBL["interior-point"])
    ax.set_xlabel("Number of pinion teeth, $z$ (-)"); ax.set_ylabel("Optimal volume (cm$^3$)")
    ax.set_xticks(z); ax.legend(frameon=False)
    save(fig, "fig3_volume_vs_z")


def fig_margins(pp):
    res = [r for r in pp["results"] if r["label"] == "optimisation"]
    names = [c["name"] for c in res[0]["constraints"]]
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    h = 0.38
    for j, r in enumerate(res):
        m = [c["margin_pct"] for c in r["constraints"]]
        ax.barh(y + (j - 0.5) * h, m, height=h - 0.04, color=C[r["method"]], label=LBL[r["method"]])
    ax.set_yticks(y); ax.set_yticklabels(names); ax.invert_yaxis()
    ax.set_xlabel("Constraint margin (% of limit), 0 = active")
    ax.axvline(0, color="k", lw=0.8)
    for j, r in enumerate(res):
        for i, c in enumerate(r["constraints"]):
            if c["active"]:
                ax.text(0.8, i + (j - 0.5) * h, "active", va="center", fontsize=8, color=C[r["method"]])
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.16), ncol=2)
    ax.grid(axis="y", visible=False)
    save(fig, "fig4_constraint_margins")


def main():
    ms = loadmat(os.path.join(RES, "multistart.mat"))
    iz = loadmat(os.path.join(RES, "integer_z.mat"))
    pp = json.load(open(os.path.join(RES, "postprocess.json")))
    post = {"best_sqp": float(np.min(ms["f_sqp"]))}
    fig_convergence(ms, post)
    fig_multistart(ms)
    fig_integer_z(iz)
    fig_margins(pp)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
