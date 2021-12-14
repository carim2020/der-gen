# Der-Gen
---

A software package to iterate over all possible substituents, generated from an already existing database. 

## Requirements
- OpenBabel 3

## Usage

1. Create a conda environment:

    `conda create -n <new_environment_name> python=3.6`
   
   `conda activate <new_environment_name>`
3. Install OpenBabel 3 with conda:

    `conda -c conda-forge install openbabel`
4. Create the following directories in the project folder:
- `in/`
- `out/`
- `out_red/`
- `error/`
4. Fill the `in/` directory with the `xyz` files from the initial database.

5. Start the script with:
    
    `python3 start.py`

6. The newly generated files are in the `out/`, `out_red/` and `error/` folders.