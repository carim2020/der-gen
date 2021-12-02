import os

from enum import Enum, auto

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

ATOMIC_SYMBOLS = dict()
with open(os.path.join(ROOT_DIR, "AtomicSymbolDict.txt"), "r") as f:
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
