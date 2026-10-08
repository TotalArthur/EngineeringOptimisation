"""solvers.py

Author: AWD Labs
Student ID: 52104479

Wrappers that run the two gradient based methods on the scaled problem.
  'sqp'            -> SciPy SLSQP        (stand in for fmincon 'sqp')
  'interior-point' -> SciPy trust-constr (stand in for fmincon 'interior-point')
Analytic-quality (complex step) gradients are supplied to both.
"""
import time
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, Bounds, BFGS
import speed_reducer_problem as prob

# tight tolerances, reported in the results (mirrors the MATLAB options)
OPTIONS = {
    "sqp": dict(ftol=1e-12, maxiter=1000, disp=False),
    "interior-point": dict(gtol=1e-10, xtol=1e-12, barrier_tol=1e-10, maxiter=3000, verbose=0),
}
FEAS_TOL = 1e-6       # max normalised violation for "feasible"
METHODS = ("sqp", "interior-point")


def free_indices(fixed):
    return [i for i in range(7) if not (fixed and i in fixed)]


def solve(method, xs0, p, fixed=None, record=False):
    """Solve from scaled start xs0 (length = number of free vars).
    Returns dict with x (display units, full), f (cm^3), viol, success, nit, nfev, time, history."""
    n = len(xs0)
    fun = lambda v: prob.objective_value(v, p, fixed)
    jac = lambda v: prob.objective_grad(v, p, fixed)
    gfun = lambda v: prob.constraint_values(v, p, fixed)
    gjac = lambda v: prob.constraint_jac(v, p, fixed)
    hist = []

    def log(v):
        if record:
            hist.append((float(fun(v)) * p["f_scale"] / p["cm"] ** 3, prob.max_violation(v, p, fixed)))

    t0 = time.perf_counter()
    if method == "sqp":
        cons = [{"type": "ineq", "fun": lambda v: -gfun(v), "jac": lambda v: -gjac(v)}]
        cb = (lambda v: log(v)) if record else None
        r = minimize(fun, xs0, jac=jac, method="SLSQP", bounds=Bounds(np.zeros(n), np.ones(n)),
                     constraints=cons, options=OPTIONS["sqp"], callback=cb)
        success = bool(r.success)
        nfev = int(r.nfev)
    elif method == "interior-point":
        nlc = NonlinearConstraint(gfun, -np.inf, 0.0, jac=gjac, hess=BFGS())
        cb = (lambda v, st: log(v) or False) if record else None
        r = minimize(fun, xs0, jac=jac, hess=BFGS(), method="trust-constr",
                     bounds=Bounds(np.zeros(n), np.ones(n)), constraints=[nlc],
                     options=OPTIONS["interior-point"], callback=cb)
        success = bool(r.success) or r.status in (1, 2)
        nfev = int(r.nfev)
    else:
        raise ValueError("solvers:badMethod unknown method '%s'" % method)
    dt = time.perf_counter() - t0
    xs = np.clip(r.x, 0.0, 1.0)
    # full display-units vector
    full = np.zeros(7)
    free = free_indices(fixed)
    xs_full = np.zeros(7)
    xs_full[free] = xs
    full = prob.from_scaled(xs_full, p)
    if fixed:
        for i, v in fixed.items():
            full[i] = v
    viol = prob.max_violation(xs, p, fixed)
    return dict(method=method, xs=xs, x=full, f=float(fun(xs)) * p["f_scale"] / p["cm"] ** 3,
                viol=viol, feasible=viol <= FEAS_TOL, success=success,
                ok=bool(success and viol <= FEAS_TOL), nit=int(r.nit), nfev=nfev,
                time=dt, history=hist, message=str(r.message))
