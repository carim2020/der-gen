import unittest

from openbabel import openbabel as ob
from typing import List, Dict
from Helper import Vector3
from Translator import dergen2ob, ob2dergen
from Definitions import ATOMIC_SYMBOLS


class TestTranslator(unittest.TestCase):
    def setUp(self) -> None:
        self.conv = ob.OBConversion()
        self.conv.SetInAndOutFormats('smi', 'svg')
        self.obm = ob.OBMol()
        self.conv.ReadString(self.obm, "C#CC1CCC1")
        self.obm.AddHydrogens()
        builder = ob.OBBuilder()
        builder.Build(self.obm)

    def test_ob2dergen(self):
        mol = ob2dergen(self.obm)

        # Check for equal number of atoms and then check all atoms
        ob_atoms: List[ob.OBAtom] = [_ for _ in ob.OBMolAtomIter(self.obm)]
        self.assertEqual(len(mol.atoms), len(ob_atoms))
        for ob_atom, ind_atom in zip(ob_atoms, mol.atoms):
            atom = mol.get_atom(ind_atom)
            self.assertEqual(Vector3(ob_atom.x(), ob_atom.y(), ob_atom.z()), atom.coord)
            self.assertEqual(ob_atom.GetAtomicNum(), atom.atomic_num)
            self.assertEqual(ob_atom.GetAtomicMass(), atom.atomic_mass)
            self.assertEqual(ATOMIC_SYMBOLS[atom.atomic_num], atom.symbol)
            self.assertEqual(ob_atom.GetIdx(), ind_atom + 1)

        # Check all of the bonds and neighbours
        # Just check the indexes. No need to create atoms as they are already revised
        ob_neighbours: Dict[int, List[int]] = {i: list() for i in range(len(ob_atoms))}
        ob_bonds: Dict[str, int] = dict()
        for b in ob.OBMolBondIter(self.obm):
            beg = b.GetBeginAtomIdx() - 1
            end = b.GetEndAtomIdx() - 1
            order = b.GetBondOrder()
            ob_bonds["{}_{}".format(beg, end)] = order
            ob_bonds["{}_{}".format(end, beg)] = order
            ob_neighbours[beg].append(end)
            ob_neighbours[end].append(beg)

        self.assertDictEqual({i: sorted([a for a in mol.neighbours[i]]) for i in mol.neighbours},
                             {i: sorted(ob_neighbours[i]) for i in ob_neighbours})
        self.assertDictEqual(ob_bonds, {key: mol.bonds[key].value for key in mol.bonds})

    def test_dergen2ob(self):
        mol = ob2dergen(self.obm)
        new_obm = dergen2ob(mol)
        conv = ob.OBConversion()
        conv.SetOutFormat("inchikey")
        old_inchi = conv.WriteString(self.obm)
        new_inchi = conv.WriteString(new_obm)
        self.assertEqual(old_inchi, new_inchi)


if __name__ == "__main__":
    unittest.main()

