# Parameters
mols_active = "mol_nums.txt"
nmol = 3375
nat_per_mol = 36



nat = nmol * at_per_mol

# Open and parse file (no error-checking)
with open(mols_active, 'r') as f:
    txt = f.read().strip("\n")
    mol_nums = txt.split()

# Create the file
s = "&ENERGY_DECOMP\n\t" + f"INDEX_MOL_DECOMP {txt}"
s += "\n\t" + f"NUM_ACTIVE_ATOMS {len(mol_nums) * at_per_mol}"
s += "\n&END ENERGY_DECOMP"
with open("DECOMP.include", 'w') as f:
    f.write(s)
