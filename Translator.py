from openbabel import openbabel as ob
from Molecule import Molecule
from Helper import Vector3
from Atom import Atom
from copy import copy
from typing import List, Dict
import os


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"
__status__ = "Prototype"


SYMBOLS = dict()
WORK_FOLDER = "."
with open(os.path.join(WORK_FOLDER, "AtomicSymbolDict.txt"), "r") as f:
    for line in f.readlines():
        (k, v) = line.split('\t')
        SYMBOLS[int(k)] = v[:-1]


def qca2ob(mol: Molecule) -> ob.OBMol:
    out = ob.OBMol()
    rekey: Dict[int, int] = dict()
    indexer = 1
    for a in mol.atoms:
        mol_atom = mol.atoms[a]
        atom = ob.OBAtom()
        atom.SetVector(mol_atom.coord.x, mol_atom.coord.y, mol_atom.coord.z)
        atom.SetAtomicNum(mol_atom.atomic_num)
        atom.SetIdx(indexer)
        rekey[mol_atom.id] = indexer
        indexer += 1
        out.AddAtom(atom)
    
    bond_set = set()
    for b in mol.bonds:
        ind = [int(i) for i in b.split('_')]
        ind_a = rekey[ind[0]]
        ind_b = rekey[ind[1]]
        bond_set.add((min(ind_a, ind_b), max(ind_a, ind_b), mol.bonds[b]))
    for b in bond_set:
        out.AddBond(b[0], b[1], b[2])
    
    return copy(out)


def ob2qca(obm: ob.OBMol) -> Molecule:
    out_mol = Molecule()
    rekey = dict()
    new_index = 0
    for atom in ob.OBMolAtomIter(obm):
        new_atom = Atom()
        new_atom.coord = Vector3(atom.x(), atom.y(), atom.z())
        new_atom.atomic_num = atom.GetAtomicNum()
        new_atom.symbol = SYMBOLS[new_atom.atomic_num]
        new_atom.atomic_mass = atom.GetAtomicMass()
        out_mol.add_atom(new_atom)
        rekey[atom.GetIdx()] = new_index
        new_index += 1
    bond_set = set()
    for bond in ob.OBMolBondIter(obm):
        beg = rekey[bond.GetBeginAtomIdx()]
        end = rekey[bond.GetEndAtomIdx()]
        order = bond.GetBondOrder()
        bond_set.add((min(beg, end), max(beg, end), order))
    for bond in bond_set:
        out_mol.add_bond(bond[0], bond[1], bond[2])
    
    return copy(out_mol)


def get_inchi_key(molecule: Molecule) -> str:

    ob_mol = qca2ob(molecule)
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("mol", "inchikey")

    res = conv.WriteString(ob_mol)
    return res


if __name__ == "__main__":
    # Testing
    conv = ob.OBConversion()
    conv.SetInAndOutFormats('smi', 'svg')

    obm = ob.OBMol()
    conv.ReadString(obm, "C1=CC=C2C=CC=CC2=C1")
    obm.AddHydrogens()
    builder = ob.OBBuilder()
    builder.Build(obm)

    mol = ob2qca(obm)
    # for atom database mol.atoms:
    #     index = mol.atoms[atom].id
    #     print("Atom's {} neighbours:".format(index))
    #     for a database mol.neighbours[index]:
    #         print(a.id)

    # obm2 = qca2ob(mol)
    # conv.WriteFile(obm2, os.path.join(WORK_FOLDER, 'out.svg'))
    print(get_inchi_key(mol))
    print("1S/C10H8/c1-2-6-10-8-4-3-7-9(10)5-1/h1-8H ")