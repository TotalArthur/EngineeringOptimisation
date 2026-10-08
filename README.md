# EG503X Assignment 1: Speed Reducer Optimisation

Author: AWD Labs
Student ID: 52104479

Minimise the total material volume (two gears plus two shafts) of a single stage speed reducer (Golinski benchmark) with gradient based methods.
Student ID 52104479 gives n = 9 and m = 7, so **T = 1500 + 100 n = 2400 Nm** and **u = 2.5 + 0.1 m = 3.2**.

## Which tool produced which results (read this first)

MATLAB with the Optimization Toolbox was **not available** in the environment used to build this folder. So:

| Part | Where | Status |
|---|---|---|
| Python port (SciPy SLSQP and trust-constr) | `code/python`, outputs in `results/python`, `figures/python` | **Run and verified. Every number in `results_summary.md` comes from here.** |
| MATLAB code (`fmincon` 'sqp' and 'interior-point') | `code/matlab`, outputs go to `results/matlab`, `figures/matlab` | Written to mirror the Python port line by line. The non-solver files (params, analysis, constraints, scaling, complex step gradients, validation) were executed in GNU Octave and reproduce the Python numbers to all printed digits. The `fmincon` driver, post processing, figure and unit test scripts have **not been executed**. Run them in MATLAB before relying on them. |

Python `trust-constr` is the closest SciPy analogue of `fmincon` interior-point, but it is not the same algorithm, so iteration counts and the tiny objective spread of the interior-point results will differ in MATLAB. Random starts also differ (`rng(1)` in MATLAB vs `numpy.random.default_rng(1)`), so individual start points are not identical across the two languages.

## Folder layout

```
code/matlab    MATLAB source
code/python    Python port (same file names, same structure)
results/python .mat and .csv from the verified runs (plus logs)
results/matlab written when you run the MATLAB driver
figures/python PNG (300 dpi) and PDF
figures/matlab written when you run make_figures.m
results_summary.md   key numbers, tables, assumptions, red flags
```

## How to run (in this order)

MATLAB (from `code/matlab`):

```matlab
validate_speed_reducer      % 1. validation table (T = 1000 and 2400 Nm)
run_optimisation            % 2. multistart, integer z study, reference case (needs Optimization Toolbox)
postprocess_design          % 3. tables, multipliers, Hessian, convexity test
make_figures                % 4. figures
results = runtests('test_speed_reducer');   % unit tests
```

Python (from `code/python`; needs numpy, scipy, matplotlib):

```bash
python3 validate_speed_reducer.py
python3 run_optimisation.py        # about 20 s, writes results/python
python3 postprocess.py
python3 literature_comparison.py   # Golinski Table 1 comparison (Python only)
python3 make_figures.py
python3 -m unittest test_speed_reducer -v
```

## Unit convention

The design vector is always in display units: `x = [b (cm), m (mm), z (-), l1 (cm), l2 (cm), d1 (cm), d2 (cm)]`.
`speed_reducer_analysis` converts to SI straight away (m, Pa, N m, m^3) and computes in SI.
Volumes are shown in cm^3, stresses in MPa and deflections in mm only for display (`a.display`).
The solver works on scaled variables `xs = (x - lb)/(ub - lb)` in [0, 1] and an objective divided by 1000 cm^3.

## Code structure

| File | Job |
|---|---|
| `speed_reducer_params` | struct of every constant, limit and bound. Cases: `optimisation`, `validation1000`, `validation2400`, `reference`. Takes an overrides struct. |
| `speed_reducer_analysis` | checked analysis (volume, stresses, deflections, geometry). Wraps `speed_reducer_core`, which also accepts complex input. |
| `speed_reducer_objective`, `speed_reducer_constraints` | scaled objective and 11 normalised constraints `c <= 0`, with complex step gradients |
| `speed_reducer_scale`, `speed_reducer_expand`, `complex_step_jac` | scaling, fixed variable handling (integer z study), derivatives |
| `solve_speed_reducer` | one solve with 'sqp' or 'interior-point' and tight tolerances |
| `validate_speed_reducer` | validation script |
| `run_optimisation` | main driver (no plotting) |
| `postprocess_design`, `make_figures` | post processing and plots, kept separate from the solver |
| `test_speed_reducer` | unit tests |

No global variables, everything passes through structs and function handles. Bad inputs raise errors with identifiers such as `speed_reducer_analysis:badSize`.

## Solver options (tight tolerances)

MATLAB `fmincon`: `OptimalityTolerance 1e-10`, `ConstraintTolerance 1e-10`, `StepTolerance 1e-12`, `MaxIterations 1000`, `MaxFunctionEvaluations 1e5`, user gradients for objective and constraints.
Python SLSQP: `ftol 1e-12`, `maxiter 1000`. Python trust-constr: `gtol 1e-10`, `xtol 1e-12`, `barrier_tol 1e-10`, `maxiter 3000`.
Both use complex step gradients (exact to machine precision) so no finite difference noise limits the tolerances.

## Changing parameters for a different problem

* Different student ID: `speed_reducer_params('optimisation', [], 'XXXXXXXX')`.
* Different loading, limits, bounds or material: edit the matching field in `speed_reducer_params` or pass an overrides struct, for example `struct('T', 3000, 'u', 3.5)`.
* Different gear geometry: change `tooth_offset`, `rim_offset`, `hub_diam_factor`, `web_fraction`, `hub_len_factor`.
* Other elastic coefficient reading: `cp_form` = `'one_minus_nu2'` (default, matches validation) or `'one_plus_nu2'` (as printed in the brief).
* New constraints: add a row to `speed_reducer_constraints` (and its name to the lists in post processing). Everything else adapts to the number of rows.
* Fixed variables (as in the integer z study): pass `fixed = struct('idx', 3, 'val', 20)` to the objective, constraints and solver.
