"""test_speed_reducer.py

Author: AWD Labs
Student ID: 52104479

Plain unittest checks: validation values, edge cases, derivative and
scaling consistency. Run:  python3 -m unittest test_speed_reducer -v
"""
import unittest
import numpy as np
from speed_reducer_params import speed_reducer_params
from speed_reducer_analysis import speed_reducer_analysis
import speed_reducer_problem as prob

X_VAL = np.array([3.5, 7.0, 22.0, 7.4, 7.8, 3.5, 5.2])


class TestValidation(unittest.TestCase):
    def setUp(self):
        self.a = speed_reducer_analysis(X_VAL, speed_reducer_params("validation2400"))["display"]

    def test_stresses_at_2400(self):
        self.assertAlmostEqual(self.a["sigma_b_MPa"], 323, delta=1)
        self.assertAlmostEqual(self.a["sigma_c_MPa"], 532, delta=1)
        self.assertAlmostEqual(self.a["sigma_s_MPa"][0], 512, delta=1)
        self.assertAlmostEqual(self.a["sigma_s_MPa"][1], 454, delta=1)

    def test_deflections_at_2400(self):
        self.assertAlmostEqual(self.a["y_mm"][0], 0.0179, delta=5e-5)
        self.assertAlmostEqual(self.a["y_mm"][1], 0.0043, delta=5e-5)

    def test_1000_does_not_match(self):
        b = speed_reducer_analysis(X_VAL, speed_reducer_params("validation1000"))["display"]
        self.assertLess(b["sigma_b_MPa"], 0.5 * 323)

    def test_volume_is_torque_independent_and_known_gap(self):
        v = self.a["V_cm3"]
        self.assertAlmostEqual(v, 3958.95, delta=0.05)     # value from the brief's formula
        self.assertGreater(abs(v - 4147) / 4147, 0.04)     # documents the unresolved 4.5 % gap


class TestEdgeCases(unittest.TestCase):
    def test_bad_size(self):
        with self.assertRaises(ValueError):
            speed_reducer_analysis(np.ones(6), speed_reducer_params())

    def test_nan(self):
        x = X_VAL.copy(); x[0] = np.nan
        with self.assertRaises(ValueError):
            speed_reducer_analysis(x, speed_reducer_params())

    def test_non_positive(self):
        x = X_VAL.copy(); x[2] = 0
        with self.assertRaises(ValueError):
            speed_reducer_analysis(x, speed_reducer_params())

    def test_bad_case(self):
        with self.assertRaises(ValueError):
            speed_reducer_params("nonsense")

    def test_student_id_parameters(self):
        p = speed_reducer_params("optimisation")
        self.assertEqual(p["T"], 2400.0)
        self.assertAlmostEqual(p["u"], 3.2)

    def test_stress_scales_with_torque(self):
        a1 = speed_reducer_analysis(X_VAL, speed_reducer_params("validation1000"))
        a2 = speed_reducer_analysis(X_VAL, speed_reducer_params("validation2400"))
        self.assertAlmostEqual(a2["sigma_b"] / a1["sigma_b"], 2.4, places=9)


class TestProblem(unittest.TestCase):
    def setUp(self):
        self.p = speed_reducer_params("optimisation")

    def test_scaling_roundtrip(self):
        xs = np.random.default_rng(0).random(7)
        np.testing.assert_allclose(prob.to_scaled(prob.from_scaled(xs, self.p), self.p), xs, atol=1e-14)

    def test_complex_step_matches_central_difference(self):
        xs = np.random.default_rng(3).random(7)
        J = prob.constraint_jac(xs, self.p)
        h = 1e-6
        Jd = np.zeros_like(J)
        for k in range(7):
            e = np.zeros(7); e[k] = h
            Jd[:, k] = (prob.constraint_values(xs + e, self.p) - prob.constraint_values(xs - e, self.p)) / (2 * h)
        np.testing.assert_allclose(J, Jd, rtol=1e-5, atol=1e-7)

    def test_known_optimum_feasible(self):
        x = np.array([3.5, 7.0, 17.0, 7.3, 7.4, 3.44361, 5.0])
        self.assertLess(prob.max_violation(prob.to_scaled(x, self.p), self.p), 1e-4)


if __name__ == "__main__":
    unittest.main()
