# TODO: Some refactoring is needed. The Molecule class has too many responsibilities.
# TODO: Make some kind of memory and speed benchmark for a substantial load (Not 3 molecules but 300)

from Atom import Atom
from Helper import Vector3
from Definitions import CycleSelection, BondType
from typing import Dict, List
from copy import copy
from openbabel import openbabel

__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


class Molecule:

    def __init__(self):
        """
            Constructor for the Molecule class.
        """
        # IMPORTANT: Everything except the dictionary for atoms uses their indexes as they are
        # unique to every atom in the molecule.
        self.__atoms: Dict[int, Atom] = dict()  # Dictionary of all the atoms' database the molecule
        self.__bonds: Dict[str, BondType] = dict()  # Dictionary of all the bonds' database the molecule
        self.__sites: List[int] = []  # Contains the atoms' database cycles
        self.__indexer = 0
        self.__neighbours: Dict[int, List[int]] = dict()  # Dictionary of lists of all the neighbours of the atoms
        self.__inchi_key: str = ""

    def add_atom(self, atom: Atom) -> int:
        self.__atoms[self.__indexer] = copy(atom)
        self.__atoms[self.__indexer].id = self.__indexer
        self.__neighbours[self.__indexer] = []
        self.__indexer += 1
        return self.__indexer - 1

    def del_atom(self, index: int) -> None:
        # Iterates over every bond in the molecule and removes the ones, connected with the
        keys = self.__bonds.keys()
        d_b = [k for k in keys if index in [int(i) for i in k.split("_")]]
        n_ind = []
        for k in keys:
            indexes = [int(i) for i in k.split('_')]
            if index in indexes:
                for ind in indexes:
                    if ind is not index:
                        n_ind.append(ind)
        for b in d_b:
            self.__bonds.pop(b)
        for ind in set(n_ind):
            self.__neighbours[ind].remove(index)
        self.__atoms.pop(index)

    # Returns reference to the atom!
    def get_atom(self, index: int) -> Atom:
        if index not in self.__atoms:
            raise KeyError("There is no atom with such id")
        return self.__atoms[index]

    def find_atom(self, coordinates: Vector3) -> Atom:
        for atom in self.__atoms:
            if self.__atoms[atom].coord == coordinates:
                return self.__atoms[atom]

    def add_bond(self, ind_a: int, ind_b: int, bond_type: BondType) -> None:
        self.__bonds["{}_{}".format(ind_a, ind_b)] = bond_type
        self.__bonds["{}_{}".format(ind_b, ind_a)] = bond_type
        self.__neighbours[ind_a].append(ind_b)
        self.__neighbours[ind_b].append(ind_a)

    def del_bond(self, ind_a: int, ind_b: int) -> None:
        self.__bonds.pop("{}_{}".format(ind_a, ind_b))
        self.__bonds.pop("{}_{}".format(ind_b, ind_a))
        self.__neighbours[ind_a].remove(ind_b)
        self.__neighbours[ind_b].remove(ind_a)

    @property
    def num_atoms(self) -> int:
        return len(self.__atoms)

    @property
    def atoms(self) -> Dict[int, Atom]:
        return copy(self.__atoms)

    @property
    def bonds(self) -> Dict[str, BondType]:
        return copy(self.__bonds)

    @property
    def sites(self) -> List[int]:
        return copy(self.__sites)

    @property
    def neighbours(self) -> Dict[int, List[int]]:
        return self.__neighbours

    @property
    def inchi_key(self) -> str:
        if self.__inchi_key == "":
            ob_mol = dergen2ob(self)
            conv = openbabel.OBConversion()
            conv.SetInAndOutFormats("mol", "inchikey")

            self.__inchi_key = conv.WriteString(ob_mol)

        return self.__inchi_key

    @inchi_key.setter
    def inchi_key(self, key: str) -> None:
        self.__inchi_key = key

    def __find_cycles(self) -> List[int]:
        passed: Dict[int, int] = {_: 0 for _ in self.__atoms}
        mark: Dict[int, bool] = {_: False for _ in self.__atoms}
        parent: Dict[int, int] = {_: -1 for _ in self.__atoms}

        if len(mark) == 0:
            raise RuntimeError("Molecule has no atoms")
        # Get the first atom
        stack = [list(self.__atoms.values())[0]]
        parent[stack[0].id] = stack[0].id
        passed[stack[0].id] = 1
        while len(stack) != 0:  # while the stack is not empty
            # get the first element of the stack
            cur = stack[-1]
            stack.pop()
            for atom_id in self.__neighbours[cur.id]:
                # Cycle check
                if passed[atom_id] == 1:
                    # print("Cycle found: {}".format(cur.id))
                    temp = cur.id
                    mark[atom_id] = True
                    mark[temp] = True
                    while temp != parent[atom_id]:
                        mark[temp] = True
                        temp = parent[temp]
                    mark[temp] = True
                    # print(self.__atoms[temp].coord)

                # Add the neighbours that are not yet discovered
                if passed[atom_id] == 0:
                    parent[atom_id] = cur.id
                    passed[atom_id] = 1
                    stack.append(self.__atoms[atom_id])

            passed[cur.id] = 2

        # Update the cycles list
        cycles: List[int] = [ind for ind, m in enumerate(mark) if m]
        return copy(cycles)

    def get_sites(self, selection: CycleSelection) -> None:
        # Why it has to be with value but not without?
        # Potential BUG
        if selection.value is CycleSelection.CYCLES.value:
            self.__sites = self.__find_cycles()

        elif selection.value is CycleSelection.NON_CYCLES.value:
            cycle = self.__find_cycles()
            self.__sites = [self.__atoms[k] for k in self.__atoms
                            if self.__atoms[k].symbol != "H" and self.__atoms[k] not in cycle]
        elif selection.value is CycleSelection.ALL.value:
            self.__sites = [self.__atoms[k] for k in self.__atoms
                            if self.__atoms[k].symbol != "H"]

    def get_neighbours(self, index: int) -> tuple:
        if index >= len(self.__sites):
            raise IndexError("Index is out of the potential sites")
        # print("Get Neighbour")
        c2 = None
        h = None
        c1 = self.__atoms[self.__sites[index]]

        for atom_id in self.__neighbours[self.__sites[index]]:
            if self.__atoms[atom_id].atomic_num != 1:
                c2 = self.__atoms[atom_id]
                # print("Carbon")
            elif self.__atoms[atom_id].atomic_num == 1:
                # print("Hydrogen")
                h = self.__atoms[atom_id]
        return c1, c2, h

    def generate_new_compound(self, c1: Atom, c2: Atom, h: Atom) -> Atom:
        """
                This method should be inherited from all the classes and used to
                form new substituent on the targeted place by changing the H atom
                for whatever substituent is desired.
        """
        raise NotImplementedError("Overwrite this function to use it with specific substituent.")

    def revert_to_original(self, c: Atom, sub: Atom) -> None:
        """
                This method reverts to the original molecule.
                Useful, will there be any need for iteration.
        """
        hydrogen = Atom()
        hydrogen.coord = sub.coord
        hydrogen.symbol = "H"
        hydrogen.atomic_num = 1

        def remove(parent: Atom, substituent: Atom):
            neigh = copy(self.neighbours[substituent.id])
            for atom_id in neigh:
                if atom_id != parent.id:
                    remove(substituent, self.__atoms[atom_id])

            if len(self.neighbours[substituent.id]) < 2:
                self.del_atom(substituent.id)

        remove(c, sub)

        new_id = self.add_atom(hydrogen)
        self.add_bond(c.id, new_id, BondType.SINGLE)

    def __eq__(self, other: 'Molecule'):
        return other.inchi_key == self.inchi_key


if __name__ == "__main__":

    import os
    from Translator import ob2dergen, dergen2ob
    from openbabel import openbabel as ob
    from Definitions import ROOT_DIR

    # Testing  code
    a = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('smi', 'svg')
    conv.ReadString(a, "C1=CC=C2C(=C1)C=CC2")
    a.AddHydrogens()
    builder = ob.OBBuilder()
    builder.Build(a)

    mol = ob2dergen(a)
    mol.get_sites(CycleSelection.CYCLES)
    print(mol.sites)
    print([mol.get_atom(key).])

    obm = dergen2ob(mol)
    conv.WriteFile(obm, os.path.join(ROOT_DIR, "original.svg"))
