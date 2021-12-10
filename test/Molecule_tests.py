
import unittest

from Atom import Atom
from Helper import Vector3
from Definitions import SiteSelection, ATOMIC_SYMBOLS, BondType
from Molecule import Molecule
from typing import List, Dict
from copy import deepcopy


class TestMolecule(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.mol = Molecule()
        self.atom_list: List[Atom] = list()
        self.bond_list: Dict[str, BondType] = dict()

        with open("../res/2-propylylcyclopentenone.xyz") as file:
            content = file.readlines()
            for line in content[2:]:
                symbol = line.split()[0]
                coord = Vector3(float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
                atomic_num = list(ATOMIC_SYMBOLS.values()).index(symbol)+1
                self.atom_list.append(Atom(coord, atomic_num, symbol))

        if not self.atom_list:
            raise FileNotFoundError("File 2-propylylcyclopentenone.xyz not found. Should be in the res/ directory.")

        self.rekey: Dict[int, int] = dict()
        for ind, a in enumerate(self.atom_list):
            self.rekey[ind] = self.mol.add_atom(a)

        for i in range(len(self.atom_list)):
            for j in range(i+1, len(self.atom_list)):
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
                elif self.atom_list[j].symbol == ATOMIC_SYMBOLS[1] and length < 1.5:
                    self.bond_list["{}_{}".format(i, j)] = BondType.SINGLE
                    self.bond_list["{}_{}".format(j, i)] = BondType.SINGLE
                    self.mol.add_bond(i, j, BondType.SINGLE)

    def takeDown(self):
        atoms = self.mol.atoms
        for a in atoms:
            self.mol.del_atom(a)
        self.atom_list.clear()
        self.bond_list.clear()
        del self.mol

    def test_init(self):
        self.assertIsInstance(self.mol, Molecule)

    def test_add_atom(self):
        empty_atom = Atom()
        index = self.mol.add_atom(empty_atom)
        self.assertEqual(index + 1, len(self.mol.atoms))
        atom_dict: Dict[int, Atom] = {self.rekey[_]: self.atom_list[_] for _ in range(len(self.atom_list))}
        atom_dict[index] = empty_atom

        self.assertListEqual(sorted([_ for _ in atom_dict]), sorted([_ for _ in self.mol.atoms]))
        for key in self.mol.atoms:
            self.assertEqual(self.mol.get_atom(key).atomic_num, atom_dict[key].atomic_num)
            self.assertTrue(self.mol.get_atom(key).coord == atom_dict[key].coord)
            self.assertEqual(self.mol.get_atom(key).symbol, atom_dict[key].symbol)

    def test_get_atom(self):
        for ind, a in enumerate(self.atom_list):
            self.assertEqual(a, self.mol.get_atom(self.rekey[ind]))

    def test_find_atom(self):
        for ind, a in enumerate(self.atom_list):
            self.assertEqual(self.rekey[ind], self.mol.find_atom(a.coord))

    def test_del_atom(self):
        ea_id = self.mol.add_atom(Atom())
        for atom_id in self.mol.atoms:
            if atom_id != ea_id:
                self.mol.add_bond(atom_id, ea_id, BondType.SINGLE)

        self.mol.del_atom(ea_id)
        self.assertNotIn(ea_id, self.mol.atoms)
        for atom_id in self.mol.atoms:
            self.assertNotIn(ea_id, self.mol.neighbours[atom_id])
            self.assertNotIn(f"{ea_id}_{atom_id}", self.mol.bonds)
            self.assertNotIn(f"{atom_id}_{ea_id}", self.mol.bonds)

        # Test whether all other atoms and bonds are intact
        atom_dict: Dict[int, Atom] = {self.rekey[_]: self.atom_list[_] for _ in range(len(self.atom_list))}
        self.assertListEqual(sorted([_ for _ in atom_dict]), sorted([_ for _ in self.mol.atoms]))
        for key in self.mol.atoms:
            self.assertEqual(self.mol.get_atom(key), atom_dict[key])

        self.assertDictEqual(self.bond_list, self.mol.bonds)

    def test_add_bond(self):
        ea_ind = self.mol.add_atom(Atom())
        self.mol.add_bond(ea_ind, self.rekey[0], BondType.SINGLE)

        self.assertIn(ea_ind, self.mol.neighbours[self.rekey[0]])
        self.assertIn(self.rekey[0], self.mol.neighbours[ea_ind])

        for b_str in self.mol.bonds:
            ind_1 = int(b_str.split("_")[0])
            ind_2 = int(b_str.split("_")[1])
            if ea_ind == ind_1:
                self.assertEqual(ind_2, self.rekey[0])
                self.assertEqual(self.mol.bonds[b_str], BondType.SINGLE)
            if ea_ind == ind_2:
                self.assertEqual(ind_1, self.rekey[0])
                self.assertEqual(self.mol.bonds[b_str], BondType.SINGLE)
            self.assertNotEqual(ind_1, ind_2)

    def test_del_bond(self):
        ea_id = self.mol.add_atom(Atom())

        for atom_id in self.mol.atoms:
            self.assertRaises(KeyError, self.mol.del_bond, ea_id, atom_id)

        for atom_id in self.mol.atoms:
            if ea_id != atom_id:
                self.mol.add_bond(ea_id, atom_id, BondType.SINGLE)

        mol_keys = list(self.mol.atoms.keys())
        for ind in mol_keys:
            atom_id = mol_keys[ind]
            if atom_id == ea_id:
                continue
            self.mol.del_bond(ea_id, atom_id)

            self.assertNotIn(atom_id, self.mol.neighbours[ea_id])
            self.assertNotIn(ea_id, self.mol.neighbours[atom_id])
            self.assertNotIn(f"{atom_id}_{ea_id}", self.mol.bonds)
            self.assertNotIn(f"{ea_id}_{atom_id}", self.mol.bonds)

            for new_id in mol_keys[ind+1:]:
                atom_id_2 = mol_keys[new_id]
                if atom_id_2 == ea_id:
                    continue
                self.assertIn(atom_id_2, self.mol.neighbours[ea_id])
                self.assertIn(ea_id, self.mol.neighbours[atom_id_2])
                self.assertIn(f"{atom_id_2}_{ea_id}", self.mol.bonds)
                self.assertIn(f"{ea_id}_{atom_id_2}", self.mol.bonds)

    def test_get_sites_ALL(self):
        sites = [self.rekey[_]
                 for _ in range(len(self.atom_list))
                 if self.atom_list[_].symbol != ATOMIC_SYMBOLS[1]]
        self.mol.get_sites(SiteSelection.ALL)
        self.assertListEqual(sorted(sites), sorted(self.mol.sites))

    def test_get_sites_CYCLES(self):
        self.mol.get_sites(SiteSelection.CYCLES)
        self.assertListEqual(sorted(self.mol.sites), [0, 1, 2])

    def test_get_sites_NON_CYCLES(self):
        self.mol.get_sites(SiteSelection.NON_CYCLES)
        self.assertListEqual(sorted(self.mol.sites), [3, 4, 5])

    def test_get_neighbours(self):
        self.mol.get_sites(SiteSelection.CYCLES)
        atoms_ids = self.mol.get_neighbours(0)
        self.assertIs(atoms_ids, None)
        c1_id, c2_id, h_id = self.mol.get_neighbours(1)

        self.assertIn(c1_id, self.mol.atoms)
        self.assertIn(c2_id, self.mol.atoms)
        self.assertIn(h_id, self.mol.atoms)

        self.assertIn(c1_id, self.mol.neighbours[c2_id])
        self.assertIn(c1_id, self.mol.neighbours[h_id])
        self.assertIn(c2_id, self.mol.neighbours[c1_id])
        self.assertIn(h_id, self.mol.neighbours[c1_id])

        self.assertIn(f"{h_id}_{c1_id}", self.mol.bonds)
        self.assertIn(f"{c1_id}_{h_id}", self.mol.bonds)
        self.assertIn(f"{c2_id}_{c1_id}", self.mol.bonds)
        self.assertIn(f"{c1_id}_{c2_id}", self.mol.bonds)

    def test_revert_to_original(self):
        self.mol.get_sites(SiteSelection.CYCLES)
        c1_id, c2_id, h_id = self.mol.get_neighbours(1)

        original_mol = deepcopy(self.mol)

        # A very large substituent with a lot of branches
        new_ids = [self.mol.add_atom(Atom(coord=self.mol.get_atom(h_id).coord)) for _ in range(10)]
        for i in range(1, 10):
            self.mol.add_bond(new_ids[i-1], new_ids[i], BondType.SINGLE)
            self.mol.add_bond(self.mol.add_atom(Atom()), new_ids[i-1], BondType.SINGLE)
            self.mol.add_bond(self.mol.add_atom(Atom()), new_ids[i-1], BondType.SINGLE)

        # New molecule
        self.mol.del_atom(h_id)
        self.mol.add_bond(c1_id, new_ids[0], BondType.SINGLE)

        self.mol.revert_to_original(c1_id, new_ids[0])
        from Translator import dergen2ob
        from openbabel import openbabel as ob
        ob_original = dergen2ob(original_mol)
        ob_mol = dergen2ob(self.mol)

        conv = ob.OBConversion()
        conv.SetOutFormat("svg")
        conv.WriteFile(ob_original, "original_mol.svg")
        conv.WriteFile(ob_mol, "mol.svg")
        self.assertEqual(original_mol, self.mol)


if __name__ == '__main__':
    unittest.main()
