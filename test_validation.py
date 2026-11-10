"""Negative controls: the replay comparison must reject the failures it exists to catch."""
import unittest
from pathlib import Path
from compare_results import text_compare, chain_compare, compare_run

HERE = Path(__file__).parent


class ValidationTests(unittest.TestCase):
    def test_changed_verdict_rejected(self):
        self.assertTrue(text_compare('variation 11.83% < 5%: FAIL', 'variation 11.83% < 5%: PASS'))

    def test_changed_value_rejected(self):
        self.assertTrue(text_compare('ln E = -1.2662', 'ln E = -1.2663'))

    def test_chain_tolerance(self):
        old = (HERE / 'checks/out/chain_checks_dr1.out').read_text()
        self.assertFalse(chain_compare(old, old))
        self.assertFalse(chain_compare(old, old.replace('w_p = -0.9452', 'w_p = -0.9453', 1)))
        self.assertTrue(chain_compare(old, old.replace('w_p = -0.9452', 'w_p = -0.9552', 1)))
        self.assertTrue(chain_compare(old, old.replace('-> gate PASS', '-> gate FAIL', 1)))

    def test_roundoff_correlations_equal(self):
        self.assertFalse(chain_compare('corr(w_p, wa) = +6.13e-16;', 'corr(w_p, wa) = -2.2e-14;'))
        self.assertTrue(chain_compare('corr(w_p, wa) = +6.13e-16;', 'corr(w_p, wa) = +3.0e-02;'))

    def test_absolute_tolerance_boundary(self):
        # The DR1 replay differs by exactly 0.001 at this printed percentile.
        for value, rejected in [('-0.7801', False), ('-0.780', False), ('-0.7799', True)]:
            for old, new in [('-0.781', value), (value, '-0.781')]:
                with self.subTest(old=old, new=new):
                    self.assertEqual(bool(chain_compare('percentile ' + old, 'percentile ' + new)), rejected)

    def test_relative_tolerance_boundary(self):
        for value, rejected in [('9999.0001', False), ('9999', False), ('9998.9999', True)]:
            for old, new in [('10000', value), (value, '10000')]:
                with self.subTest(old=old, new=new):
                    self.assertEqual(bool(chain_compare('count ' + old, 'count ' + new)), rejected)

    def test_deterministic_tolerance_boundary(self):
        for value, rejected in [('9.999e-10', False), ('1e-9', False), ('1.0001e-9', True)]:
            with self.subTest(value=value):
                self.assertEqual(bool(text_compare('value 0', 'value ' + value)), rejected)

    def test_published_dr1_percentile_boundary(self):
        run = {'out': 'out/chain_checks_dr1.out'}
        old = (HERE / 'checks' / run['out']).read_text()
        self.assertIn('99.9th weighted percentile -0.781;', old)
        for value, rejected in [('-0.780', False), ('-0.779', True)]:
            with self.subTest(value=value):
                new = old.replace('99.9th weighted percentile -0.781;',
                                  '99.9th weighted percentile ' + value + ';', 1)
                self.assertEqual(bool(compare_run(run, HERE / '_no_such_work_copy', new)), rejected)

    def test_reference_is_the_published_output(self):
        # End to end through compare_run: the reference must come from the package, not from a work copy
        # the check may have rewritten, so a small numeric change in stdout is still caught.
        run = {'out': 'out/prior_volume.out'}
        old = (HERE / 'checks' / run['out']).read_text()
        bogus = HERE / '_no_such_work_copy'
        self.assertFalse(compare_run(run, bogus, old))
        self.assertTrue(compare_run(run, bogus, old.replace('0.281887', '0.281888', 1)))


if __name__ == '__main__':
    unittest.main()
