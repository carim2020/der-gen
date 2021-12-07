import math
import os

from Definitions import ATOMIC_SYMBOLS, CycleSelection, BondType
from Molecule import Molecule
from Atom import Atom
from Helper import cross_product


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


class Propylyl(Molecule):
    def __init__(self, molecule: Molecule = Molecule()):
        super().__init__()
        for a_i in molecule.atoms:
            self.add_atom(molecule.atoms[a_i])

        bond_set = set()
        for b_i in molecule.bonds:
            index = [int(i) for i in b_i.split('_')]
            bond_set.add((min(index[0], index[1]), max(index[0], index[1]), molecule.bonds[b_i]))

        for b in bond_set:
            self.add_bond(b[0], b[1], b[2])

    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int):
        """
                Extension of the Molecule class to create propylyl compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """

        bond_vector = self.get_atom(h_id).coord - self.get_atom(c1_id).coord
        bond_vector = bond_vector / bond_vector.length()

        norm_vector = cross_product(self.get_atom(h_id).coord - self.get_atom(c1_id).coord,
                                    self.get_atom(c2_id).coord - self.get_atom(c1_id).coord)
        norm_vector = norm_vector / norm_vector.length()

        new_c1 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c1.coord = self.get_atom(c1_id).coord + bond_vector * 1.5  # a bit shorter than sp3 bond
        new_c1_id = self.num_atoms

        new_c2 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c2.coord = new_c1.coord + bond_vector * 1.37  # sp-sp bond
        new_c2_id = self.num_atoms + 1

        new_c3 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c3.coord = new_c2.coord + bond_vector * 1.46  # sp-sp3 bond bond
        new_c3_id = self.num_atoms + 2

        new_H_vec = bond_vector * math.cos(math.radians(70)) + norm_vector * math.sin(math.radians(70))

        new_H = [Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1])]
        for ind, hydro in enumerate(new_H):
            hydro.coord = new_c3.coord + new_H_vec.rotate_around_axis(bond_vector, 2/3*math.pi * ind)

        self.del_atom(h_id)
        self.add_atom(new_c1)
        self.add_atom(new_c2)
        self.add_atom(new_c3)
        new_H_id = [self.add_atom(_) for _ in new_H]

        self.add_bond(c1_id, new_c1_id, BondType.SINGLE)
        self.add_bond(new_c1_id, new_c2_id, BondType.TRIPLE)
        self.add_bond(new_c2_id, new_c3_id, BondType.SINGLE)
        for ind in new_H_id:
            self.add_bond(new_c3.id, ind, BondType.SINGLE)

        return new_c1_id


if __name__ == "__main__":

    from os.path import join
    from Translator import ob2dergen, dergen2ob
    from openbabel import openbabel as ob

    WORK_FOLDER = os.path.join(os.getcwd(), 'out')

    # Testing code
    a = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('smi', 'svg')
    conv.ReadString(a, "C1=CC=C2C(=C1)C=CC2")
    a.AddHydrogens()
    builder = ob.OBBuilder()
    builder.Build(a)

    mol = Propylyl(ob2dergen(a))
    mol.get_sites(CycleSelection.CYCLES)
    # print([(m.id, m.symbol) for m in mol.cycles])
    c1, c2, h = mol.get_neighbours(0)

    cn = mol.generate_new_compound(c1, c2, h)

    obm = dergen2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, 'out_2.svg'))

    mol.revert_to_original(c1, cn)
    obm = dergen2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, "original.svg"))
