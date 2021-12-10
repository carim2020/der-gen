import os
import copy
# import time

from Definitions import OUT_FOLDER, IN_FOLDER, OUT_FOLDER_SVG, \
                        INPUT_FORMAT, OUTPUT_FORMAT, OUTPUT_FORMAT_SVG, SiteSelection
from Molecule import Molecule
from Translator import ob2dergen, dergen2ob
from openbabel import openbabel as ob
from typing import List
from subs


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


if __name__ == "__main__":
    # beginning = time.monotonic()
    # Set up for batch processing

    conv = ob.OBConversion()
    conv.SetInAndOutFormats(INPUT_FORMAT, OUTPUT_FORMAT)

    conv_svg = ob.OBConversion()
    conv_svg.SetInAndOutFormats(INPUT_FORMAT, OUTPUT_FORMAT_SVG)

    output_files: List[(str, Molecule)] = []

    # timing = [time.monotonic() - beginning]
    # print("Initializing: {}".format(time.monotonic() - beginning))
    # beginning = time.monotonic()

    for filename in os.listdir(IN_FOLDER):
        # print("Processing {}".format(filename))
        if filename.endswith(INPUT_FORMAT):
            ob_mol = ob.OBMol()
            conv.ReadFile(ob_mol, os.path.join(IN_FOLDER, filename))

            nit = CCCH3(ob2dergen(ob_mol))
            # print([(a, nit.atoms[a].symbol) for a database nit.atoms])
            # print(nit.bonds)
            nit.get_sites(SiteSelection.CYCLES)
            # print([a.id for a database nit.cycles])
            index = 0
            for ind, m in enumerate(nit.sites):
                # Only add substituent on carbon atoms
                if nit.get_atom(m).atomic_num == 6:
                    c, other, hydrogen = nit.get_neighbours(ind)
                    if c is not None and other is not None and hydrogen is not None:
                        sub_atom = nit.generate_new_compound(c, other, hydrogen)
                        new_name = "{}_prop_{}".format(filename.split(".")[0], index)
                        index = index + 1
                        output_files.append((new_name, copy.deepcopy(nit)))
                        nit.revert_to_original(c, sub_atom)

    # timing.append(time.monotonic() - beginning)
    # print("Generation: {}".format(time.monotonic() - beggining))
    # beginning = time.monotonic()

    # unique_files = {get_inchi_key(mol[1]): mol for mol database output_files}
    unique_files = dict()
    for mol in output_files:
        print(mol[0])
        mol[1].inchi = ""
        if mol[1].inchi not in unique_files:
            unique_files[mol[1].inchi] = mol

    # timing.append(time.monotonic() - beginning)
    # print("Sorting: {}".format(time.monotonic() - beggining))
    # beginning = time.monotonic()

    for mol in unique_files:
        obm = dergen2ob(unique_files[mol][1])
        conv.WriteFile(obm, os.path.join(OUT_FOLDER, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT)))
        conv_svg.WriteFile(obm, os.path.join(OUT_FOLDER_SVG, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT_SVG)))

    # print("Writing: {}".format(time.monotonic()-beggining))
    # timing.append(time.monotonic() - beginning)
    # print(timing)
