from Reduced import Reduced
from Atom import Atom
from Definitions import BondType, ATOMIC_SYMBOLS


class NNN(Reduced):
    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int) -> int:
        """
                Extension of the Molecule class to create nitrile compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """
        unit_v = (self.get_atom(h_id).coord - self.get_atom(c1_id).coord).unit()

        newN1 = Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7)
        newN1.coord = self.get_atom(c1_id).coord + unit_v * 1.25

        newN2 = Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7)
        newN2.coord = newN1.coord + unit_v * 1.25

        newN3 = Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7)
        newN3.coord = newN2.coord + unit_v * 1.25

        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h_id)
        newN1_id = self.add_atom(newN1)
        newN2_id = self.add_atom(newN2)
        newO_id = self.add_atom(newN3)

        self.add_bond(c1_id, newN1_id, BondType.SINGLE)
        self.add_bond(newN1_id, newN2_id, BondType.DOUBLE)
        self.add_bond(newN2_id, newO_id, BondType.DOUBLE)

        # self.end_modify()

        return newN1_id
