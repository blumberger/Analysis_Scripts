import coord_wrap

from src.wrappers import wrap_utils


def calc_RDF(pos, ats_per_mol, atom_types, ABC):
	"""
	A wrapper for the calc RDF C function.

	This is supposed to make the C function a bit easier to use so one
	doesn't need to construct the input dict etc...

	Inputs:
		* pos <list> => The atomic positions in shape (natom, 3)
		* ats_per_mol <int> => How many atoms in a molecule
		* atom_types <list> => The atomic types in shape (natom)
		* ABC <list> => The cell vectors:
		         ABC = [ [self.Var['xlo'], self.Var['xhi'], self.Var['xy']],
        		         [self.Var['ylo'], self.Var['yhi'], self.Var['xz']],
                		 [self.Var['zlo'], self.Var['zhi'], self.Var['yz']] ]
	Outputs:
		<arr> wrapped coords.
	"""
    # Handle the atom type parameters
	unique_types = np.unique(atom_types)
	at_type_dict = {T: i for i, T in enumerate(unique_types)}
	for T in at_type_dict:
		atom_types[atom_types == T] = at_type_dict[T]
	atom_types = atom_types.astype(int)

	# Get the molecular indices
	nmol = wrap_utils.get_nmol(len(pos), ats_per_mol)
	mol_inds, _ = np.mgrid[0:nmol, 0:ats_per_mol]
	mol_inds = mol_inds.flatten()

	# Create standard python objects out of everything
	pos = wrap_utils.to_list(pos)
	mol_inds = wrap_utils.to_list(mol_inds)
	atom_types = wrap_utils.to_list(atom_types)

	coord_wrap.calc_RDF({'pos': pos, 'mol_nums': mol_inds, 'atom_types': atom_types,
				 			   'ABC': ABC})

	# return np.array(wrapped_coords)