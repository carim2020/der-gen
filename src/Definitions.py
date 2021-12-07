import os

from enum import Enum, auto

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

ATOMIC_SYMBOLS = dict()
with open(os.path.join(ROOT_DIR, "../res/AtomicSymbolDict.txt"), "r") as f:
    for line in f.readlines():
        (k, v) = line.split('\t')
        ATOMIC_SYMBOLS[int(k)] = v[:-1]

OUT_FOLDER = os.path.join(ROOT_DIR, "out")
OUT_FOLDER_SVG = os.path.join(ROOT_DIR, "out_svg")
IN_FOLDER = os.path.join(ROOT_DIR, "in")
OUTPUT_FORMAT = "xyz"
OUTPUT_FORMAT_SVG = "svg"
INPUT_FORMAT = "xyz"


class CycleSelection(Enum):
    CYCLES = auto()
    NON_CYCLES = auto()
    ALL = auto()


class BondType(Enum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 5

    @staticmethod
    def from_ob_bond_bype(bond_type: int) -> 'BondType':
        if bond_type == BondType.SINGLE.value:
            return BondType.SINGLE
        elif bond_type == BondType.DOUBLE.value:
            return BondType.DOUBLE
        elif bond_type == BondType.TRIPLE.value:
            return BondType.TRIPLE
        elif bond_type == BondType.AROMATIC.value:
            return BondType.AROMATIC
