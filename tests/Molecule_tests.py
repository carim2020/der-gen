import unittest

from Atom import Atom
from Helper import Vector3
from Definitions import CycleSelection, ATOMIC_SYMBOLS, BondType
from Molecule import Molecule
from typing import List, Dict


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.mol = Molecule()
        self.atom_list: List[Atom] = list()
        self.bond_list: Dict[str, BondType] = dict()

        with open("test.xyz") as file:
            content = file.readlines()
            for line in content[2:]:
                symbol = line.split()[0]
                coord = Vector3(float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
                atomic_num = list(ATOMIC_SYMBOLS.values()).index(symbol)
                self.atom_list.append(Atom(coord, atomic_num, symbol))

        if not self.atom_list:
            raise FileNotFoundError("File test.xyz not found. Should be in the test/ directory.")

        for ind, a in enumerate(self.atom_list):
            self.atom_list[ind].id = self.mol.add_atom(a)

        for i in range(len(self.atom_list)):
            for j in range(i, len(self.atom_list)):
                length = (self.atom_list[i].coord - self.atom_list[j].coord).length()
                if 1.5 < length < 2:
                    self.bond_list["{}_{}".format(i, j)] = BondType.SINGLE
                    self.bond_list["{}_{}".format(j, i)] = BondType.SINGLE
                    self.mol.add_bond(i, j, BondType.SINGLE)
                elif 1.4 < length < 1.5:
                    self.bond_list["{}_{}".format(i, j)] = BondType.DOUBLE
                    self.bond_list["{}_{}".format(j, i)] = BondType.DOUBLE
                    self.mol.add_bond(i, j, BondType.DOUBLE)
                elif 1.1 < length < 1.3:
                    self.bond_list["{}_{}".format(i, j)] = BondType.TRIPLE
                    self.bond_list["{}_{}".format(j, i)] = BondType.TRIPLE
                    self.mol.add_bond(i, j, BondType.TRIPLE)

    def takeDown(self):
        atoms = self.mol.atoms
        for a in atoms:
            self.mol.del_atom(a)
        self.atom_list.clear()
        self.bond_list.clear()
        del self.mol

    def test_init(self):
        self.assertIsInstance(self.mol, Molecule)

    def test_add(self):
        empty_atom = Atom()
        index = self.mol.add_atom(empty_atom)
        self.assertEqual(index + 1, len(self.mol.atoms))
        atom_dict: Dict[int, Atom] = {a.id: a for a in self.atom_list}
        atom_dict[index] = empty_atom

        self.assertListEqual(sorted([_ for _ in atom_dict]), sorted([_ for _ in self.mol.atoms]))
        for key in self.mol.atoms:
            self.assertEqual(self.mol.get_atom(key).atomic_num, atom_dict[key].atomic_num)
            self.assertTrue(self.mol.get_atom(key).coord == atom_dict[key].coord)
            self.assertEqual(self.mol.get_atom(key).symbol, atom_dict[key].symbol)


if __name__ == '__main__':
    unittest.main()
