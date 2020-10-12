# Parameters
AOM_COEFF_mol_file = "./src/data/AOM_COEFFs/Pentacene_Lammps_Single_Mol.include"
AOM_COEFF_sing_at = "./src/data/AOM_COEFFs/single_inert_at.include"

mols_active = "mol_nums.txt"

nmol = 3375
at_per_mol = 36



# Open and parse file (no error-checking)
with open(mols_active, 'r') as f:
    ltxt = f.read()
    mol_nums = tuple(map(int, ltxt.split()))

# Create the mol mask
mask = [False] * nmol
for i in mol_nums:  mask[i] = True

# Create mol strings for inert and active mols
with open(AOM_COEFF_sing_at, 'r') as f:
    txt = f.read()
    inert_mol = txt * at_per_mol

with open(AOM_COEFF_mol_file, 'r') as f:
    act_mol = f.read()

# Create the AOM_COEFF file (simple loop)
filetxt = ""
for i in mask:
    if i: filetxt += act_mol
    else: filetxt += inert_mol

with open("AOM_COEFF.include", 'w') as f:
    f.write(filetxt)

