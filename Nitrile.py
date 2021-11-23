import os

from Molecule import Molecule, CycleSelection
from Atom import Atom


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"
__status__ = "Prototype"


class Nitrile(Molecule):
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
                Extension of the Molecule class to create nitrile compounds
                @param c1: Origin
                @param c2: Child of C1
                @param h: Child of C1
        """
        newC = Atom()
        newC.atomic_num = 6
        newC.symbol = "C"
        # newC.id = self.num_atoms
        newC.coord = h.coord
        
        newN = Atom()
        newN.atomic_num = 7
        newN.symbol = "N"
        # newN.id = self.num_atoms + 1
        newN.coord = (h.coord * 2) - c1.coord
        
        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h.id)
        newC.id = self.add_atom(newC)
        newN.id = self.add_atom(newN)
        print([a.id for a in self.neighbours[c1.id]])
        
        self.add_bond(c1.id, newC.id, 1)
        self.add_bond(newC.id, newN.id, 3)
        print([a.id for a in self.neighbours[c1.id]])
        print([a.id for a in self.neighbours[newC.id]])

        # self.end_modify()

        return self.atoms[newC.id]


if __name__ == "__main__":

    from os.path import join
    from Translator import ob2qca, qca2ob
    from openbabel import openbabel as ob

    WORK_FOLDER = os.getcwd()

    # Testing code
    a = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('smi', 'svg')
    conv.ReadString(a, "C1=CC=C2C(=C1)C=CC2")
    a.AddHydrogens()
    builder = ob.OBBuilder()
    builder.Build(a)

    mol = Nitrile(ob2qca(a))
    mol.get_sites(CycleSelection.CYCLES)
    c1, c2, h = mol.get_neighbours(0)

    cn = mol.generate_new_compound(c1, c2, h)

    obm = qca2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, 'out_2.svg'))

    mol.revert_to_original(c1, cn)
    obm = qca2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, "original.svg"))
