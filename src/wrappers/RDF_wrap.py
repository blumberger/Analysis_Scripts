import numpy as np

# The RDF C module
import RDF

from src.wrappers import wrap_utils


def calc_RDF(pos, ats_per_mol, atom_types, ABC,
			 types_to_calc_1='all', types_to_calc_2='all', cutoff=False, dr=False):
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
		OPTIONAL:
		* types_to_calc_1 <list | 'all'> => Which type to calculate the RDF for (first in pair)
		* types_to_calc_2 <list | 'all'> => Which type to calculate the RDF for (second in pair)
		* cutoff <float> => A cutoff for the RDF calculation (doesn't affect runtime)
		* dr <float> => The spacing between bins (doesn't affect runtime)

	Outputs:
		<arr>, <arr> radii and rdf
	"""
	# Deal with optional arguments
	x_len = np.max(pos[:, 0]) - np.min(pos[:, 0])
	y_len = np.max(pos[:, 1]) - np.min(pos[:, 1])
	z_len = np.max(pos[:, 2]) - np.min(pos[:, 2])
	if cutoff is False:	cutoff = min([x_len, y_len, z_len]) / 2.
	if dr is False:	dr = cutoff / 1000

    # Handle the atom type parameters
	unique_types = np.unique(atom_types)
	at_type_dict = {T: i for i, T in enumerate(unique_types)}
	for T in at_type_dict:
		atom_types[atom_types == T] = at_type_dict[T]
	atom_types = atom_types.astype(int)

	# Handle which atoms to calculate for
	if types_to_calc_1 == 'all': types_to_calc_1 = list(unique_types[:])
	if types_to_calc_2 == 'all': types_to_calc_2 = list(unique_types[:])
	for i, T in enumerate(types_to_calc_1): 	types_to_calc_1[i] = at_type_dict[T]
	for i, T in enumerate(types_to_calc_2):		types_to_calc_2[i] = at_type_dict[T]

	# Get the molecular indices
	nmol = wrap_utils.get_nmol(len(pos), ats_per_mol)
	mol_inds, _ = np.mgrid[0:nmol, 0:ats_per_mol]
	mol_inds = mol_inds.flatten()

	# Create standard python objects out of everything
	pos = wrap_utils.to_list(pos)
	mol_inds = wrap_utils.to_list(mol_inds)
	atom_types = wrap_utils.to_list(atom_types)

	rdf, radii = RDF.calc_RDF({'pos': pos, 'mol_nums': mol_inds, 'atom_types': atom_types,
				 			   'cutoff': cutoff, 'dr': dr, 'type_list1': types_to_calc_1,
				 			   'type_list2': types_to_calc_2, 'ABC': ABC})

	return np.array(radii), np.array(rdf)
