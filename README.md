# EG503X Assignment 1: Speed Reducer Optimisation

Author: AWD Labs
Student ID: 52104479 (so T = 2400 Nm and u = 3.2)

## How to run

1. Open `speed_reducer_main.m` in MATLAB (needs the Optimization Toolbox).
2. Press **Run**. It takes a few minutes.

That is all. The one file does everything, in numbered sections:

| Section | What it does |
|---|---|
| 0 | settings and every parameter, limit and bound (change things here) |
| 1 | validation against the brief's validation case (T = 1000 and 2400 Nm) |
| 2 | one solve with fmincon `sqp` and `interior-point` |
| 3 | multistart, 100 random starts per method (`rng(1)`) |
| 4 | integer teeth study, z = 17 to 28 |
| 5 | post processing: constraint and bound tables, multipliers, Hessian, convexity test, sanity check |
| 6 | reference case (u = 3, T = 1000 Nm) and comparison with Golinski (1970) |
| 7 | figures |

The functions it uses (`analysis`, `objective`, `constraints` and a few helpers) are at the bottom of the same file.
Everything is printed in the Command Window and logged to `results/run_log.txt`.
Tables and results are saved to a `results` folder, and figures (PNG at 300 dpi and PDF) to a `figures` folder, both created next to the script.

## Units

The design vector is always `x = [b (cm); m (mm); z; l1 (cm); l2 (cm); d1 (cm); d2 (cm)]`.
The `analysis` function converts to SI (m, Pa, N m, m^3) straight away and works in SI.
Results are shown in cm^3, MPa and mm.

## Changing the problem

Everything lives in section 0 at the top:

* different student ID: change `studentID` (this sets T and u),
* different loading, material, limits or bounds: edit the matching line of `p`,
* different gear geometry: edit the `dst_off`, `dw_off`, `dp_fac`, `web` and `hub_len` lines,
* more or fewer random starts: `nStarts`, `nStartsZ`, `nStartsRef`,
* new constraint: add a row in `constraints()` (and a name in `conNames`, and a row in `physical_constraints()`).

## Solver settings

`fmincon` with central finite difference gradients, `OptimalityTolerance 1e-8`, `ConstraintTolerance 1e-8`, `StepTolerance 1e-10`, `MaxIterations 1000`. They are set in `run_solver`.

## Other files

* `results_summary.md`: key numbers, tables, assumptions and red flags for the report.
* `python_check/`: optional. A Python version of the same calculation that was used to cross-check the numbers (MATLAB was not available when the code was written). You do not need it. Run `python3 run_optimisation.py` there if you want to regenerate its results.

## How the MATLAB script was tested

MATLAB itself was not available. The full script was run end to end in GNU Octave, with `fmincon` temporarily replaced by Octave's own `sqp` solver, and it reproduced the Python results (volume 3018.17 cm^3, same multipliers, same Hessian eigenvalues). So the logic, tables and file writing are tested, but the real `fmincon` calls have not been run. If MATLAB shows a small error on first run, that is the most likely place. Differences to expect with the real `fmincon`: the interior-point result will sit very slightly inside the feasible region (a few thousandths of a cm^3 above SQP), and iteration counts and run times will differ.
