# Results summary: speed reducer optimisation

Author: AWD Labs
Student ID: 52104479 (n = 9, m = 7), so T = 2400 Nm and u = 3.2

**Provenance.** All numbers below come from the Python port in `code/python` (SciPy SLSQP as the SQP method and SciPy trust-constr as the interior-point method), which was actually executed. MATLAB `fmincon` was not available, so the MATLAB files in `code/matlab` are unexecuted for the solver parts (see README). Octave did run the MATLAB analysis, constraint and validation code and matched Python exactly. Times are for the machine used here and are not MATLAB times.

## 1. Validation

Inputs: u = 3, q = 2.54, kv = 2.1, E = 200 GPa, nu = 0.3, z = 22, m = 7.0 mm, b = 3.5 cm, l1 = 7.4 cm, l2 = 7.8 cm, d1 = 3.5 cm, d2 = 5.2 cm. Contact stress uses cp = E / (2 pi (1 - nu^2)), see 1.2.

| Quantity | Target | T = 1000 Nm | Error | T = 2400 Nm | Error |
|---|---|---|---|---|---|
| Volume (cm^3) | 4147 | 3959 | -4.53 % | 3959 | -4.53 % |
| sigma_b (MPa) | 323 | 134.6 | -58.32 % | 323.1 | +0.04 % |
| sigma_c (MPa) | 532 | 343.5 | -35.43 % | 532.2 | +0.03 % |
| sigma_s1 (MPa) | 512 | 213.5 | -58.30 % | 512.4 | +0.08 % |
| sigma_s2 (MPa) | 454 | 189.1 | -58.35 % | 453.8 | -0.03 % |
| y1 (mm) | 0.0179 | 0.007442 | -58.42 % | 0.01786 | -0.22 % |
| y2 (mm) | 0.0043 | 0.001789 | -58.40 % | 0.004293 | -0.17 % |

(The -0.22 % and -0.17 % deflection differences are just rounding of the 4 figure targets: 0.01786 rounds to 0.0179 and 0.004293 rounds to 0.0043.)

### 1.1 Which torque reproduces the targets

**T = 2400 Nm reproduces every stress and deflection target. T = 1000 Nm does not**: every torque dependent quantity comes out 58 % low (a factor of 1/2.4). The brief's "T = 1000 Nm" for the validation case therefore looks like a typo, or the validation numbers were generated with T = 2400 Nm. 2400 Nm is also exactly what the optimisation formula gives for a student ID ending in 9 (1500 + 100 x 9), which suggests the validation sheet was produced with that ID. I have not changed any constant to get this match.

### 1.2 Contact stress formula (sigma_c)

The PDF shows sigma_c = (cp * kv T / (b m^2 z^2) * (1+u)/u)^(1/2) with z squared in the denominator and cp = E / (2 pi (1 + nu^2)). Tested at T = 2400 Nm:

| Reading of cp | sigma_c (MPa) | Error vs 532 |
|---|---|---|
| E / (2 pi (1 - nu^2)) | 532.15 | +0.03 % |
| E / (2 pi (1 + nu^2)) (as printed) | 486.23 | -8.60 % |

Only the 1 - nu^2 reading matches, and it is also the physically standard Hertz form, so the printed "1 + nu^2" is most likely a typo. The 1 - nu^2 form is used throughout, with a switch (`cp_form`) to flip it. Using the printed form instead does not change the optimum (sigma_c is not active, 624.3 MPa vs 683.3 MPa, both below 800 MPa).

### 1.3 Volume does NOT match (unresolved, 4.5 %)

The brief's volume formula (implemented exactly as printed) gives 3959 cm^3 against the 4147 cm^3 target. The volume does not depend on torque, so this is not the T issue. What I checked:

* The formula is the standard Golinski benchmark objective: it gives the same 3958.9 cm^3 (to 0.1) as the usual closed form 0.7854 b m^2 (3.3333 z^2 + 14.9334 z - 43.0934) - 1.508 b (d1^2 + d2^2) + 7.4777 (d1^3 + d2^3) + 0.7854 (l1 d1^2 + l2 d2^2), and about 2994.4 cm^3 at the commonly quoted benchmark optimum design (3.5, 0.7, 17, 7.3, 7.7153, 3.3505, 5.2867). That design is from my memory of the benchmark, not from the course PDF.
* No single input has a sensible value that explains 4147 (matching would need u = 3.13, z = 22.85, m = 7.25 mm, b = 3.78 cm, and so on).
* A brute force search over variants of the geometry constants (offsets, web fraction, hub length, shaft length) found only multi-parameter coincidences with no physical meaning. I did not adopt any.

The brief's formula is the authoritative one, so it is used. All volumes in this report are on that basis. **This should be raised with the course team.** The stresses and deflections are not affected.

## 2. Optimisation case

T = 2400 Nm, u = 3.2, q = 2.54, kv = 2.1, E = 200 GPa, nu = 0.3, limits sigma_b 650, sigma_c 800, sigma_s 550 MPa, y 0.075 mm, bounds as in the brief.

Solver options (tight tolerances). SQP (SLSQP): ftol 1e-12, maxiter 1000. Interior point (trust-constr): gtol 1e-10, xtol 1e-12, barrier_tol 1e-10, maxiter 3000. Complex step gradients for both. Feasible means max normalised constraint violation <= 1e-6. A run counts as successful if the solver reports success and the result is feasible. MATLAB equivalents are listed in the README.

### 2.1 Best design

| Variable | SQP | Interior point | Unit | Bound status |
|---|---|---|---|---|
| b | 3.500 | 3.500 | cm | free (pinned by b/m = 5) |
| m | 7.000 | 7.000 | mm | lower bound |
| z | 17.00 | 17.00 | - | lower bound |
| l1 | 7.300 | 7.300 | cm | lower bound |
| l2 | 7.400 | 7.400 | cm | free (pinned by 1.1 d2 + 1.9) |
| d1 | 3.444 | 3.444 | cm | free (pinned by shaft 1 stress) |
| d2 | 5.000 | 5.000 | cm | lower bound |
| **Volume** | **3018 cm^3** (3018.174) | **3018 cm^3** (3018.176) | | |

The interior-point result sits about 1e-5 inside the feasible region, as expected from a barrier method, which explains its 0.002 cm^3 higher volume.

### 2.2 Constraints at the SQP optimum (active tolerance: relative physical margin < 1e-4)

| Constraint | Value | Limit | Margin | Status | Multiplier* |
|---|---|---|---|---|---|
| Gear bending stress (MPa) | 418.2 | 650 | 35.66 % | inactive | 0 |
| Gear contact stress (MPa) | 683.3 | 800 | 14.59 % | inactive | 0 |
| Shaft 1 combined stress (MPa) | 550.0 | 550 | 0 % | **active** | 0.3089 |
| Shaft 2 combined stress (MPa) | 545.4 | 550 | 0.840 % | inactive (close) | 0 |
| Shaft 1 deflection (mm) | 0.02368 | 0.075 | 68.43 % | inactive | 0 |
| Shaft 2 deflection (mm) | 0.005550 | 0.075 | 92.60 % | inactive | 0 |
| Shaft 1 min length 1.5 d1 + 1.9 (cm) | 7.065 | 7.300 (= l1) | 3.214 % | inactive | 0 |
| Shaft 2 min length 1.1 d2 + 1.9 (cm) | 7.400 | 7.400 (= l2) | 0 % | **active** | 0.1433 |
| b/m lower limit (5) | 5.000 | b/m | 0 % | **active** | 1.163 |
| b/m upper limit (12) | 5.000 | 12 | 58.33 % | inactive | 0 |
| Overall size m z (1+u) (cm) | 49.98 | 160 | 68.76 % | inactive | 0 |

*Multipliers are for the scaled problem (variables in [0,1], objective divided by 1000 cm^3, constraints normalised), so only their sign and relative size are meaningful. In the Python port they are computed from the KKT conditions by non-negative least squares on the active set, not read from the solver. The MATLAB post processing uses the `fmincon` multipliers directly. The interior-point design gives the same active set and multipliers to 5 digits. Margin for the b/m row is quoted against the actual b/m ratio, and likewise the length rows against l1 and l2.

Bound multipliers (scaled): m 0.7213, z 2.097, l1 0.01403, d2 0.2939, all positive.

### 2.3 Optimality and structure

* Three active constraints plus four active bounds make 7 active, linearly independent conditions for 7 variables, so the optimum is a **vertex**. All 7 multipliers are strictly positive, so the first-order (KKT) conditions hold with strict complementarity and the point is a strict local minimum. The null space is empty, so the second-order condition is satisfied trivially.
* First-order optimality measure (max norm of the Lagrangian gradient, scaled space): 2.2e-16 (SQP), 4.4e-16 (interior point), which is just the least squares residual with the active set fixed. It is not an independent solver output.
* Hessian of the Lagrangian (scaled space, central differences of complex step gradients) eigenvalues: -0.6414, -0.1281, -0.0002710, 0.0005710, 0.05686, 0.5669, 2.105. The matrix is **indefinite** in the full space. That does not contradict the local result because the reduced Hessian on the empty null space has no eigenvalues.
* Sanity check: all 11 original, unnormalised constraints and all bounds are satisfied at both best designs (checked from the physical quantities, independently of the normalised rows).

### 2.4 Convexity assessment

Midpoint convexity test along 5000 random segments inside the bounds (violation means h(midpoint) > average of the end values):

| Function | Segments violating | Comment |
|---|---|---|
| Objective (volume) | 1727 (34.5 %) | clearly **non-convex** |
| Overall size m z (1+u) | 2552 (51.0 %) | non-convex (product of m and z) |
| Shaft 2 deflection | 83 (1.7 %) | slightly non-convex |
| Shaft 1 deflection | 46 (0.9 %) | slightly non-convex |
| Shaft 1 / 2 combined stress | 2 / 1 (<0.1 %) | almost convex, tiny gaps |
| Gear bending, gear contact stress | 0 | no violation found |
| Four linear geometric rows | 0 (gaps ~1e-14, numerical noise) | convex (they are linear) |

**What the evidence supports:** the problem is non-convex (objective and overall size constraint fail the test, Lagrangian Hessian is indefinite), so no global optimum can be guaranteed by theory. A convexity test by sampling can only prove non-convexity, never prove convexity.
**What it does not support:** it does not show there is more than one local minimum. In practice all 200 multistart runs ended at the same design (section 3), and the optimum is a nondegenerate vertex. That is strong empirical evidence for global optimality, not proof.

## 3. Multistart (100 random starts per method, rng seed 1, uniform in the bounds, z continuous)

| Statistic | SQP | Interior point |
|---|---|---|
| Success rate (flag success and feasible) | 100 / 100 | 100 / 100 |
| Best objective (cm^3) | 3018.174 | 3018.176 |
| Worst objective (cm^3) | 3018.174 | 3018.399 |
| Spread worst - best (cm^3) | 3.1e-10 | 0.2232 |
| Mean (cm^3) / std | 3018.174 / 4.7e-11 | 3018.192 / 0.0322 |
| Iterations, mean (median) | 9.59 (10) | 27.58 (28) |
| Function evaluations, mean | 9.66 | 20.88 |
| Time per run, mean (total) | 0.0199 s (1.99 s) | 0.1261 s (12.61 s) |

Every successful run of both methods lands on the same vertex (the same active set). The interior-point spread of 0.223 cm^3 (0.007 %) is the barrier method stopping slightly inside the feasible region, not different local optima. Only 12 of the 100 interior-point values are within 1e-6 relative of their own best for the same reason, so I judge agreement by design variables and active set rather than by that count.

Observation: in an earlier run with a different (variable) normalisation of the four linear geometric constraints, 4 of 100 SQP runs stopped with "positive directional derivative for linesearch" at the same optimum (violation below 6e-10). After switching those rows to constant scaling, all 100 succeeded. So SQP success can depend on constraint scaling at these tight tolerances.

## 4. Integer teeth

Fixed z = 17..28, re-optimised the other six variables (20 random starts per z and method; all feasible except 1 of 20 SQP starts at z = 28, which did not affect the best value).

| z | Optimal volume (cm^3), SQP | Interior point |
|---|---|---|
| 17 | **3018** | 3018 |
| 18 | 3214 | 3214 |
| 19 | 3420 | 3420 |
| 20 | 3637 | 3637 |
| 21 | 3864 | 3864 |
| 22 | 4101 | 4101 |
| 23 | 4348 | 4348 |
| 24 | 4606 | 4606 |
| 25 | 4873 | 4873 |
| 26 | 5151 | 5151 |
| 27 | 5439 | 5439 |
| 28 | 5737 | 5737 |

Volume rises monotonically with z. The continuous optimum already sits at z = 17, the lower bound and an integer, so the **best integer design is the continuous optimum and the rounding penalty is 0 cm^3 (0 %)**. If z were forced to 18 instead, the penalty would be 195.9 cm^3 (6.49 %). This is a feature of these bounds and this loading, not a general result.

## 5. Literature comparison (Golinski 1970, Table 1)

Table 1 values were checked against the scanned paper (f 2236.35 and 2247.79 cm^3, designs as in the brief; module is given in cm in the paper). His loads, material limits and bounds differ from the course case, so only trends are compared. Reference case = u = 3, T = 1000 Nm with the course limits and bounds (50 starts per method, all successful, SQP best 2697.154 cm^3, same vertex for both methods).

| Design | f reported / this model (cm^3) | b | m (mm) | z | l1 | l2 | d1 | d2 | Notes |
|---|---|---|---|---|---|---|---|---|---|
| Golinski crude Monte Carlo | 2236.35 / 2671.9 | 4.4 | 6.0 | 17 | 7.3 | 8.1 | 3.4 | 5.0 | m below our bound |
| Golinski stray process | 2247.79 / 3049.9 | 3.6 | 7.0 | 18 | 6.6 | 8.2 | 2.8 | 5.2 | l1 below our bound, b/m = 5.14 |
| Reference case (T = 1000, u = 3) | 2697 (this model) | 3.500 | 7.000 | 17 | 7.300 | 7.400 | 2.800 | 5.000 | active: l2 length, b/m = 5 |
| **This work (T = 2400, u = 3.2)** | **3018** (this model) | 3.500 | 7.000 | 17 | 7.300 | 7.400 | 3.444 | 5.000 | active: shaft 1 stress, l2 length, b/m = 5 |

Trends that agree: smallest allowed number of teeth (17 or 18), smallest allowed module, small pinion shaft diameter, l1 and d2 at or near their lower limits, and b/m at or near its lower limit of 5 (stray process 5.14). In our case b ends up at 3.5 cm only because b/m = 5 with m = 7 mm.
Effect of the larger torque and ratio: going from (T = 1000, u = 3) to (T = 2400, u = 3.2) raises the volume by 321.0 cm^3 (11.9 % above the reference case). The only design change is d1 growing from its lower bound 2.800 to 3.444 cm, because shaft 1 stress becomes active.

**Red flag:** the objective as printed in the brief, evaluated at Golinski's two Table 1 designs, gives 2672 and 3050 cm^3, not 2236 and 2248 cm^3. I cannot explain that with the scan quality of the paper's formula (the hub and shaft terms are hard to read), and I did not try to force it. It may be that his published f used slightly different geometry, so the numeric comparison of f is not meaningful, only the trends above.

## 6. Assumptions, ambiguities and things I was unsure about

1. **Validation torque.** T = 2400 Nm reproduces the targets, T = 1000 Nm does not (section 1.1). Treated as a typo in the brief.
2. **Elastic coefficient.** Used 1 - nu^2 instead of the printed 1 + nu^2, because only it matches 532 MPa (section 1.2). The optimum is unaffected.
3. **Volume validation fails by 4.5 %** and is unresolved (section 1.3). The brief's formula is used unchanged.
4. **z squared in the contact stress.** The PDF text extraction was ambiguous, but the page image shows z^2 in the denominator and that matches the target.
5. **Units of b/m.** Taken as dimensionless with b and m in the same unit (b = 3.5 cm and m = 7 mm gives exactly 5, matching the validation case sitting on the limit). Overall size m z (1+u) <= 160 uses cm.
6. **Constraint normalisation.** Stress, deflection and overall size rows are value/limit - 1. The four linear geometric rows use a constant reference (smallest allowed length or width) instead of dividing by a variable, so they stay linear and convex. This slightly departs from "g/g_max - 1" for those four rows only.
7. **Python is not MATLAB.** SLSQP and trust-constr stand in for `fmincon` sqp and interior-point. Multipliers in the Python port come from a KKT least squares on the active set. Random starts differ between languages.
8. **Active tolerance.** A constraint is called active when its relative physical margin is below 1e-4 and a variable is at a bound when within 1e-4 of its range. This is looser than the solver tolerance because the interior-point result stops about 1e-5 inside the feasible set.
9. **Convexity test** is sampling based (5000 segments, seed 2, midpoint test, relative tolerance 1e-9). It can show non-convexity but cannot prove convexity.
10. **Shaft 2 stress margin is only 0.84 %**, so it is nearly active. A small change in loading could make it active.
11. **Multistart is evidence, not proof,** of global optimality (non-convex problem). Only z fixed integer values 17 to 28 and continuous z were studied.
12. **Unexecuted MATLAB.** `run_optimisation.m`, `solve_speed_reducer.m`, `postprocess_design.m`, `make_figures.m` and `test_speed_reducer.m` were written carefully but not run. Expect to fix small syntax issues on first run, and expect interior-point results from `fmincon` to differ slightly from trust-constr.
13. Golinski's own f values could not be reproduced from the brief's formula (section 5).

## 7. Notes for the report, mapped to the assignment brief (EG503X/Y Assignment 1)

The brief marks four things: formulation (20 %), coding (30 %), solution and analysis (30 %), presentation (20 %). The report is capped at 10 pages including everything except the title page, with numbered equations and defined variables with units. This section lists what the results already cover and what still needs writing up. The realism comments are my engineering judgement, not outputs of the code.

### Task 1, formulation: points to state in the report

* Objective: total material volume of two gears plus two shafts. Design variables: b, m, z, l1, l2, d1, d2 (units in section 2.1). Constraints: 11 inequalities (gear bending and contact stress, two shaft stresses, two deflections, two minimum shaft lengths, two b/m limits, overall size) plus 7 bounds. Parameters: T, u, q, kv, E, nu and the four limits.
* Simplifications and reformulations to mention: z treated as continuous first and integer afterwards; gear 2 has u z teeth, which is not an integer in general (54.4 teeth here); variables and objective scaled to [0, 1] and 1000 cm^3 for conditioning; the four linear geometric constraints normalised by a constant to stay linear; complex step gradients instead of finite differences; no safety factors, one load case, constant q, kv, E and nu; cp taken with 1 - nu^2 (section 1.2).

### Task 2, coding: what is covered

Parameters, analysis, objective and constraints are separate files, nothing is hard coded elsewhere, there are no globals, and inputs are checked. Validation is in section 1, including the two unresolved discrepancies (T typo, volume 4.5 % off), which the report should state openly rather than hide. Unit tests cover the validation values and edge cases. The MATLAB solver scripts still need one run in MATLAB (README).

### Task 3, solution: points to state in the report

The brief asks for two methods, a comparison with the literature, convexity, active constraints, bound status and realism of the optimum.

* **Two methods and literature:** sections 2 to 5. The brief says "input parameters from the literature", while the problem sheet fixes T and u from the student ID. The report should say the main case uses the sheet's ID-based values and that the u = 3, T = 1000 Nm reference case is the closer analogue of Golinski's setting. Even that case uses different limits and bounds from his, so only trends are compared.
* **Convex or not, active constraints:** sections 2.2 to 2.4.
* **Variables at limits:** m, z, l1 and d2 sit at their lower bounds, and b, l2 and d1 are fixed by active constraints (b/m = 5, 1.1 d2 + 1.9, shaft 1 stress). So no variable is truly free. The solution is determined by bounds and constraints, which is why every start finds it. The largest bound multipliers (scaled) are for z (2.097) and m (0.7213), so relaxing those two bounds would pay off most. This is a statement about the sheet's bounds, not a general design rule.
* **Are the values realistic (judgement):**
  * z = 17 is the usual minimum pinion tooth count for a 20 degree full depth gear before undercutting, so it is realistic but right at the practical limit.
  * m = 7 mm is a standard (second choice) module, but it sits on the lower bound.
  * The 54.4 teeth on gear 2 (u z = 3.2 x 17) cannot be built. Real gears would need an integer pair such as 17 and 54 or 55, which changes the ratio slightly (3.176 or 3.235).
  * Shaft 1 runs at its full 550 MPa allowable with no safety factor, and b/m = 5 is the sheet's minimum rather than a design preference.
  * Volume is the only objective. Cost, fatigue, bearings, housing and manufacturability are not modelled.
* **What this says about the problem and practical value:** it is a small, bound dominated, vertex-type problem. Gradient methods solve it easily and reliably, and the real value of the result is in the active set and multipliers (which limits drive the design), more than in the single volume figure.
