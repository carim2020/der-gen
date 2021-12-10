from Reduced import Reduced
from Atom import Atom
from Definitions import BondType, ATOMIC_SYMBOLS


class SCN(Reduced):
    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int) -> int:
        """
                Extension of the Molecule class to create nitrile compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """
        unit_v = (self.get_atom(h_id).coord - self.get_atom(c1_id).coord).unit()

        newS = Atom(symbol=ATOMIC_SYMBOLS[16], atomic_number=16)
        newS.coord = self.get_atom(c1_id).coord + unit_v * 1.80

        newC = Atom(symbol=ATOMIC_SYMBOLS[6], atomic_number=6)
        newC.coord = newS.coord + unit_v * 1.80
        
        newN = Atom(symbol=ATOMIC_SYMBOLS[7], atomic_number=7)
        newN.coord = newC.coord + unit_v * 1.2
        
        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h_id)
        newS_id = self.add_atom(newS)
        newC_id = self.add_atom(newC)
        newN_id = self.add_atom(newN)

        self.add_bond(c1_id, newS_id, BondType.SINGLE)
        self.add_bond(newC_id, newS_id, BondType.SINGLE)
        self.add_bond(newC_id, newN_id, BondType.TRIPLE)

        # self.end_modify()

        return newS_id
