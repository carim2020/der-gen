import unittest

from subs.SCN import SCN
from Atom import Atom
from Definitions import ATOMIC_SYMBOLS, BondType
from subs_test.subs_test_template import TestSubTemplate
from Helper import Vector3


class TestThiocyanide(TestSubTemplate):
    def test_init(self):
        red = SCN(self.mol)
        self.assertIsInstance(red, SCN)

    def test_generate_new_compound(self):
        red = SCN(self.mol)
        red.generate_new_compound(1, 0, 6)
        
        self.assertEqual(len(red.atoms), 10)
        self.assertEqual(red.get_atom(8), Atom(symbol=ATOMIC_SYMBOLS[16], atomic_number=16,
                                               coord=red.get_atom(1).coord + self.unit_v * 1.8))
        self.assertEqual(red.get_atom(9), Atom(symbol=ATOMIC_SYMBOLS[6], atomic_number=6,
                                               coord=red.get_atom(8).coord + self.unit_v * 1.8))
        self.assertEqual(red.get_atom(10), Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7,
                                                coord=red.get_atom(9).coord + self.unit_v * 1.47))
        
        self.assertEqual(red.bonds["1_8"], BondType.SINGLE)
        self.assertEqual(red.bonds["8_1"], BondType.SINGLE)
        self.assertEqual(red.bonds["8_9"], BondType.SINGLE)
        self.assertEqual(red.bonds["9_8"], BondType.SINGLE)
        self.assertEqual(red.bonds["9_10"], BondType.TRIPLE)
        self.assertEqual(red.bonds["10_9"], BondType.TRIPLE)


if __name__ == "__main__":
    unittest.main()
