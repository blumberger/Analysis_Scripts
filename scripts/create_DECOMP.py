# Parameters
mols_active = "NVE_folder/Restraint_Stuff/slice_mol_nums.txt"
#mols_active = "all"
nmol = 0
at_per_mol = 36


import os


# Open and parse file (no error-checking)
if os.path.isfile(mols_active):
    with open(mols_active, 'r') as f:
        txt = f.read().strip("\n")
        mol_nums = txt.split()
else:
    nat = nmol * at_per_mol
    mol_nums = tuple(range(nmol))

txt = ' '.join(map(lambda x: str(int(x)+1), mol_nums))

# Create the file
s = "&ENERGY_DECOMP\n\t" + f"INDEX_MOL_DECOMP {txt}"
s += "\n\t" + f"NUM_ACTIVE_ATOMS {len(mol_nums) * at_per_mol}"
s += "\n&END ENERGY_DECOMP"
with open("/home/matt/Data/Work/NewPentaceneSlab/24_25_5/Fixed_Atoms_FSSH/DECOMP.include", 'w') as f:
    f.write(s)
