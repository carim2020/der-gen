import math

from Definitions import ATOMIC_SYMBOLS, BondType
from Reduced import Reduced
from Atom import Atom
from Helper import cross_product


class CCCH3(Reduced):
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

        new_c1 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c1.coord = self.get_atom(c1_id).coord + bond_vector * 1.5  # a bit shorter than sp3 bond

        new_c2 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c2.coord = new_c1.coord + bond_vector * 1.37

        new_c3 = Atom(atomic_number=6, symbol=ATOMIC_SYMBOLS[6])
        new_c3.coord = new_c2.coord + bond_vector * 1.46  # sp-sp3 bond bond

        new_H_vec = bond_vector * math.cos(math.radians(70)) + norm_vector * math.sin(math.radians(70))

        new_H = [Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1]),
                 Atom(atomic_number=1, symbol=ATOMIC_SYMBOLS[1])]
        for ind, hydro in enumerate(new_H):
            hydro.coord = new_c3.coord + new_H_vec.rotate_around_axis(bond_vector, 2/3*math.pi * ind)

        self.del_atom(h_id)
        new_c1_id = self.add_atom(new_c1)
        new_c2_id = self.add_atom(new_c2)
        new_c3_id = self.add_atom(new_c3)
        new_H_id = [self.add_atom(_) for _ in new_H]

        self.add_bond(c1_id, new_c1_id, BondType.SINGLE)
        self.add_bond(new_c1_id, new_c2_id, BondType.TRIPLE)
        self.add_bond(new_c2_id, new_c3_id, BondType.SINGLE)
        for ind in new_H_id:
            self.add_bond(new_c3_id, ind, BondType.SINGLE)

        return new_c1_id
