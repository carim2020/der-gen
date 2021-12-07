import unittest

from Reduced import Reduced
from Translator import ob2dergen
from Atom import Atom
from Definitions import CycleSelection
from openbabel import openbabel as ob
from typing import Dict


class TestReduced(unittest.TestCase):
    def setUp(self) -> None:
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        obm = ob.OBMol()
        conv.ReadString(obm, "C1=CC(=O)C(=CC1=O)CO")  # Methyl-p-benzoquinone
        obm.AddHydrogens()
        builder = ob.OBBuilder()
        builder.Build(obm)
        self.mol = ob2dergen(obm)

    def test_init(self):
        red = Reduced(self.mol)
        self.assertIsInstance(red, Reduced)

    def test_reduce(self):
        red = Reduced(self.mol)

        bl1 = red.get_atom(3).coord - red.get_atom(2).coord
        bl1 = bl1 / bl1.length() * 2.
        li_1 = Atom(coord=red.get_atom(3).coord + bl1, symbol="Li", atomic_number=3)
        li_1.id = len(red.atoms)

        bl2 = red.get_atom(7).coord - red.get_atom(6).coord
        bl2 = bl2 / bl2.length() * 2
        li_2 = Atom(coord=red.get_atom(7).coord + bl2, symbol="Li", atomic_number=3)
        li_2_id = len(red.atoms) + 1

        old_atoms: Dict[int, Atom] = red.atoms
        old_atoms[li_1.id] = li_1
        old_atoms[li_2.id] = li_2

        red.get_sites(CycleSelection.CYCLES)
        red.reduce()

        self.assertListEqual(sorted([_ for _ in old_atoms]), sorted([_ for _ in red.atoms]))
        for key in red.atoms:
            self.assertEqual(red.atoms[key].atomic_num, old_atoms[key].atomic_num)
            self.assertTrue(red.atoms[key].coord == old_atoms[key].coord)
            self.assertEqual(red.atoms[key].symbol, old_atoms[key].symbol)

if __name__ == "__main__":
    unittest.main()
