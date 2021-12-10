import unittest

from subs.CN import CN
from subs_test.subs_test_template import TestSubTemplate
from Atom import Atom
from Helper import Vector3
from Definitions import ATOMIC_SYMBOLS, BondType


class TestNitrile(TestSubTemplate):
    def test_init(self):
        red = CN(self.mol)
        self.assertIsInstance(red, CN)

    def test_generate_new_compound(self):
        red = CN(self.mol)
        red.generate_new_compound(1, 0, 6)

        self.assertEqual(len(red.atoms), 9)
        self.assertEqual(red.get_atom(8), Atom(symbol=ATOMIC_SYMBOLS[6], atomic_number=6,
                                               coord=Vector3(-1.9884643835538665, -0.8483596826630728, 0.9720840138877235)))
        self.assertEqual(red.get_atom(9), Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7,
                                               coord=Vector3(-2.873585852121078, -0.567382946071825, 1.633304143175962)))

        self.assertEqual(red.bonds["8_9"], BondType.TRIPLE)
        self.assertEqual(red.bonds["1_8"], BondType.SINGLE)


if __name__ == "__main__":
    unittest.main()
