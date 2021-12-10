import unittest

from subs.CCCH3 import CCCH3
from subs_test.subs_test_template import TestSubTemplate


class TestPropyl(TestSubTemplate):
    def test_init(self):
        red = CCCH3(self.mol)
        self.assertIsInstance(red, CCCH3)

    def test_generate_new_compound(self):
        red = CCCH3(self.mol)
        red.generate_new_compound(1, 0, 6)

        raise NotImplementedError


if __name__ == "__main__":
    unittest.main()
