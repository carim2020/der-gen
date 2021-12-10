import os

from Reduced import Reduced
from Atom import Atom
from Definitions import BondType


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


class CN(Reduced):

    def generate_new_compound(self, c1_id: int, c2_id: int, h_id: int) -> int:
        """
                Extension of the Reduced class to create nitrile compounds
                @param c1_id: Origin
                @param c2_id: Child of C1
                @param h_id: Child of C1
        """
        unit_v = (self.get_atom(h_id).coord - self.get_atom(c1_id).coord).unit()
        newC = Atom(symbol="C", atomic_number=6)
        newC.coord = self.get_atom(c1_id).coord + unit_v * 1.36
        
        newN = Atom(symbol="N", atomic_number=7)
        newN.coord = newC.coord + unit_v * 1.14
        
        # First called BeginModify to stop reindexing at every step
        # self.start_modify()

        self.del_atom(h_id)
        newC_id = self.add_atom(newC)
        newN_id = self.add_atom(newN)

        self.add_bond(c1_id, newC_id, BondType.SINGLE)
        self.add_bond(newC_id, newN_id, BondType.TRIPLE)

        # self.end_modify()

        return newC_id
