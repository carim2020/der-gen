
from Molecule import Molecule
from Atom import Atom
from Helper import Vector3
from Definitions import BondType, ATOMIC_SYMBOLS


class Reduced(Molecule):
    def __init__(self, molecule: Molecule):
        super().__init__()
        for a_i in molecule.atoms:
            self.add_atom(molecule.atoms[a_i])

        bond_set = set()
        for b_i in molecule.bonds:
            index = [int(i) for i in b_i.split('_')]
            bond_set.add((min(index[0], index[1]), max(index[0], index[1]), molecule.bonds[b_i]))

        for b in bond_set:
            self.add_bond(b[0], b[1], b[2])

    def generate_new_compound(self, c1: Atom, c2: Atom, h: Atom) -> Atom:
        raise NotImplementedError

    def reduce(self, atomic_number: int = 3, bond_length: float = 2.) -> int:
        sites = self.sites
        if not sites:
            raise RuntimeError("There are no sites in this molecule. First generate them.")

        counter: int = 0
        for atom_id in sites:
            if len(self.neighbours[atom_id]) == 3:
                for neigh_id in self.neighbours[atom_id]:
                    if self.get_atom(neigh_id).atomic_num == 8 and len(self.neighbours[neigh_id]) == 1:

                        li = Atom(symbol=ATOMIC_SYMBOLS[atomic_number], atomic_number=atomic_number)
                        bond_v: Vector3 = self.get_atom(neigh_id).coord - self.get_atom(atom_id).coord
                        bond_v = bond_v / bond_v.length() * bond_length
                        li.coord = bond_v + self.get_atom(neigh_id).coord
                        li.id = self.add_atom(li)
                        self.add_bond(neigh_id, li.id, BondType.SINGLE)

                        counter += 1

        return counter

    def oxidise(self, atomic_number: int = 3) -> None:
        atoms = self.atoms.keys()
        for a in atoms:
            if self.get_atom(a).atomic_num == atomic_number:
                self.del_atom(a)
