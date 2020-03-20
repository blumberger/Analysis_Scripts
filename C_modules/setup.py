from distutils.core import setup, Extension

RDF_module = Extension('RDF',
                        sources = ['src/RDF.c'])

rdf_description = """Python package to calculate the radial distribution function.

To use the utility use the function RDF.calc_RDF() as an argument supply a dictionary
with the following keys:
	* pos:        The positions. These should be given as a 2D list of shape (natom, 3).
	* ABC:        The vector that define the simulation box\nThe input should be of the form:\n[    [xlo, xhi, xy],\n     [ylo, yhi, xz],\n     [zlo, zhi, yz] ]
	* mol_nums    The molecular number of each atom should be 1D list of length (natom, )
	* atom_types  The type of each atom, should be 1D list of length (natom, )
	* cutoff      The cutoff distance used when calculating the RDF (reducing won't speed up the code)
	* dr          The spacing between bins in the histogram (increasing won't speed up the code)
	* type_list1  The atom types to calculate the RDF for (first index of atomic pair) list of arbitrary length
	* type_list2  The atom types to calculate the RDF for (second index of atomic pair) list of arbitrary length

"""

setup(name = "RDF",
      version = "1.0",
      description = rdf_description,
      ext_modules = [RDF_module],
      url="https://en.wikipedia.org/wiki/Radial_distribution_function",
      author="Matt Ellis",
      author_email="95ellismle@gmail.com")
