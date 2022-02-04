import os
import copy
# import time

from Definitions import OUT_FOLDER, IN_FOLDER, RED_OUT_FOLDER, \
                        INPUT_FORMAT, OUTPUT_FORMAT, OUTPUT_FORMAT_SVG, SiteSelection
from Translator import ob2dergen, dergen2ob
from openbabel import openbabel as ob
from typing import List, Tuple
from Reduced import Reduced


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


if __name__ == "__main__":
    output_files: List[Tuple[str, Reduced]] = []
    reduced_files: List[Tuple[str, Reduced]] = []
    errors: List[Tuple[str, Reduced]] = []

    conv = ob.OBConversion()
    conv.SetInFormat(INPUT_FORMAT)

    # These lines load all the modules from the subs start.pydirectory dynamically
    import subs
    from subs import *
    for name in subs.__all__:
        # print(locals()[name])
        for filename in os.listdir(IN_FOLDER):
            if filename.endswith(INPUT_FORMAT):
                ob_mol = ob.OBMol()
                conv.ReadFile(ob_mol, os.path.join(IN_FOLDER, filename))

                mol = locals()[name](ob2dergen(ob_mol))

                index = 0
                mol.get_sites(SiteSelection.ALL)
                for ind, m in enumerate(mol.sites):
                    # Only add substituent on carbon atoms
                    if mol.get_atom(m).atomic_num == 6:
                        for possible_site in mol.get_neighbours(ind):
                            if not any(_ == -1 for _ in possible_site):
                                sub_atom = mol.generate_new_compound(possible_site[0], possible_site[1], possible_site[2])
                                new_name = "{}_{}_{}".format(filename.split(".")[0], name.split(".")[-1], index)

                                output_files.append((new_name, copy.deepcopy(mol)))
                                if mol.reduce() == 2:
                                    reduced_files.append((f"Li_{new_name}", copy.deepcopy(mol)))
                                else:
                                    errors.append((f"Li_{new_name}", copy.deepcopy(mol)))
                                mol.oxidise()
                                mol.revert_to_original(possible_site[0], sub_atom)
                                index = index + 1

    unique_files = dict()
    for mol in output_files:
        mol[1].inchi_key = ""
        if mol[1].inchi_key not in unique_files:
            unique_files[mol[1].inchi_key] = mol

    for mol in unique_files:
        obm = dergen2ob(unique_files[mol][1])
        conv.WriteFile(obm, os.path.join(OUT_FOLDER, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT)))
        # conv.WriteFile(obm, os.path.join(OUT_FOLDER_SVG, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT_SVG)))

    unique_files = dict()
    for mol in reduced_files:
        mol[1].inchi_key = ""
        if mol[1].inchi_key not in unique_files:
            unique_files[mol[1].inchi_key] = mol

    for mol in unique_files:
        obm = dergen2ob(unique_files[mol][1])
        conv.WriteFile(obm, os.path.join(RED_OUT_FOLDER, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT)))
        # conv.WriteFile(obm, os.path.join(RED_OUT_FOLDER_SVG, "{}.{}".format(unique_files[mol][0], OUTPUT_FORMAT_SVG)))

    for mol in errors:
        obm = dergen2ob(mol[1])
        conv.WriteFile(obm, os.path.join(OUT_FOLDER, "{}.{}".format(mol[0], OUTPUT_FORMAT)))
        # conv.WriteFile(obm, os.path.join(OUT_FOLDER_SVG, "{}.{}".format(mol[0], OUTPUT_FORMAT_SVG)))

