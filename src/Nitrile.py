import os

from Molecule import Molecule, CycleSelection
from Atom import Atom
from Definitions import BondType


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


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

    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int) -> int:
        """
                Extension of the Molecule class to create nitrile compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """
        newC = Atom(symbol="C", atomic_number=6)
        # newC.id = self.num_atoms
        newC.coord = self.get_atom(h_id).coord
        
        newN = Atom()
        newN.atomic_num = 7
        newN.symbol = "N"
        # newN.id = self.num_atoms + 1
        newN.coord = (self.get_atom(h_id).coord * 2) - self.get_atom(c1_id).coord
        
        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h_id)
        newC_id = self.add_atom(newC)
        newN_id = self.add_atom(newN)

        self.add_bond(c1_id, newC_id, BondType.SINGLE)
        self.add_bond(newC_id, newN_id, BondType.TRIPLE)

        # self.end_modify()

        return newC_id


if __name__ == "__main__":

    from os.path import join
    from Translator import ob2dergen, dergen2ob
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

    mol = Nitrile(ob2dergen(a))
    mol.get_sites(CycleSelection.CYCLES)
    c1, c2, h = mol.get_neighbours(0)

    cn = mol.generate_new_compound(c1, c2, h)

    obm = dergen2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, 'out_2.svg'))

    mol.revert_to_original(c1, cn)
    obm = dergen2ob(mol)
    conv.WriteFile(obm, os.path.join(WORK_FOLDER, "original.svg"))
