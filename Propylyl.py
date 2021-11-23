import math
import os

from Molecule import Molecule, CycleSelection
from Atom import Atom
from Helper import cross_product, dot_product, Vector3


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"
__status__ = "Prototype"


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

    def generate_new_compound(self, c1: Atom, c2: Atom, h: Atom):
        """
                Extension of the Molecule class to create propylyl compounds
                @param c1: Origin
                @param c2: Child of C1
                @param h: Child of C1
        """

        bond_vector = h.coord - c1.coord
        bond_vector = bond_vector / bond_vector.length()

        norm_vector = cross_product(h.coord - c1.coord, c2.coord - c1.coord)
        norm_vector = norm_vector / norm_vector.length()

        new_c1 = Atom()
        new_c1.atomic_num = 6
        new_c1.symbol = "C"
        new_c1.id = self.num_atoms
        new_c1.coord = c1.coord + bond_vector * 1.5  # a little bit shorter than sp3 bond
        
        new_c2 = Atom()
        new_c2.atomic_num = 6
        new_c2.symbol = "C"
        new_c2.id = self.num_atoms + 1
        new_c2.coord = new_c1.coord + bond_vector * 1.37  # sp-sp bond

        new_c3 = Atom()
        new_c3.atomic_num = 6
        new_c3.symbol = "C"
        new_c3.id = self.num_atoms + 2
        new_c3.coord = new_c2.coord + bond_vector * 1.46  # sp-sp3 bond bond

        new_H_vec = bond_vector * math.cos(math.radians(70)) + norm_vector * math.sin(math.radians(70))

        new_H = [Atom(), Atom(), Atom()]
        for ind, hydro in enumerate(new_H):
            hydro.atomic_num = 1
            hydro.symbol = "H"
            hydro.id = self.num_atoms + ind + 3
            hydro.coord = new_c3.coord + new_H_vec.rotate_around_axis(bond_vector, 2/3*math.pi * ind)

        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h.id)
        self.add_atom(new_c1)
        self.add_atom(new_c2)
        self.add_atom(new_c3)
        for i in new_H:
            self.add_atom(i)

        self.add_bond(c1.id, new_c1.id, 1)
        self.add_bond(new_c1.id, new_c2.id, 3)
        self.add_bond(new_c2.id, new_c3.id, 1)
        for i in new_H:
            self.add_bond(new_c3.id, i.id, 1)
        
        # self.end_modify()

        return self.find_atom(new_c1.coord)


if __name__ == "__main__":

    from os.path import join
    from Translator import ob2qca, qca2ob
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

    mol = Propylyl(ob2qca(a))
    mol.get_sites(CycleSelection.CYCLES)
    # print([(m.id, m.symbol) for m in mol.cycles])
    c1, c2, h = mol.get_neighbours(0)

    cn = mol.generate_new_compound(c1, c2, h)

    obm = qca2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, 'out_2.svg'))

    mol.revert_to_original(c1, cn)
    obm = qca2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, "original.svg"))
