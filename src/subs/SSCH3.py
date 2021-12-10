import math

from Definitions import ATOMIC_SYMBOLS, BondType
from Reduced import Reduced
from Atom import Atom
from Helper import cross_product


class SSCH3(Reduced):
    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int):
        """
                Extension of the Reduced class to create propylyl compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """

        bond_vector = (self.get_atom(h_id).coord - self.get_atom(c1_id).coord).unit()

        norm_vector = cross_product(self.get_atom(h_id).coord - self.get_atom(c1_id).coord,
                                    self.get_atom(c2_id).coord - self.get_atom(c1_id).coord).unit()

        newS1 = Atom(atomic_number=16, symbol=ATOMIC_SYMBOLS[16])
        newS1.coord = self.get_atom(c1_id).coord + bond_vector * 1.8

        newS2 = Atom(atomic_number=16, symbol=ATOMIC_SYMBOLS[6])
        newS2.coord = newS1.coord + bond_vector * 1.8  # sp-sp bond

        newC3 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        newC3.coord = newS2.coord + bond_vector * 1.8  # sp-sp3 bond bond

        new_H_vec = bond_vector * math.cos(math.radians(70)) + norm_vector * math.sin(math.radians(70))

        new_H = [Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1])]
        for ind, hydro in enumerate(new_H):
            hydro.coord = newC3.coord + new_H_vec.rotate_around_axis(bond_vector, 2/3*math.pi * ind)

        self.del_atom(h_id)
        newS1_id = self.add_atom(newS1)
        newS2_id = self.add_atom(newS2)
        newC3_id = self.add_atom(newC3)
        new_H_id = [self.add_atom(_) for _ in new_H]

        self.add_bond(c1_id, newS1_id, BondType.SINGLE)
        self.add_bond(newS1_id, newS2_id, BondType.SINGLE)
        self.add_bond(newS2_id, newC3_id, BondType.SINGLE)
        for ind in new_H_id:
            self.add_bond(newC3_id, ind, BondType.SINGLE)

        return newS1_id
