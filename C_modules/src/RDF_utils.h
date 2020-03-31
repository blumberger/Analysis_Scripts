#include <Python.h>
#include <math.h>
#include <string.h>
#define cmax_length 1000

#define pi 3.14159265

/*
    A type to hold position file_data
*/
static struct Pos {
	double *x, *y, *z;
};

/*
	A type to hold 1 position.
*/
static struct Single_Pos {
	double x, y, z;
};

/*
 	A convienient way to hold the cell vectors
*/
static struct Matrix {
	struct Single_Pos a;
	struct Single_Pos b;
	struct Single_Pos c;
};

/*
    A type to hold all the file_data from a lammps file.
*/
struct Lammps_Dump {
    struct Pos R;
    struct Pos sR;
    int *mol, *type;
    int nx, ny, nz;
    int natom;
    int timestep;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound;
    int is_triclinic;
    double xy, xz, yz;
    int exit_code;
};

/*
	A struct to store the data from each type section in from the config file
*/
static struct Config_File_Type_Section {
	int N, *list;
};

/*
	A type to hold the config file data
*/
struct Config_File {
	struct Config_File_Type_Section type_1;
	struct Config_File_Type_Section type_2;
	int exit_code;
};

/*
	A type to contain which atoms belong to which lists.
*/
static struct belongs_to_lists {
	int *list1, *list2;
    int N1, N2;
};

/*
Will hold the parameters required to calculate the RDF
*/
struct Params {
	double dr;
	double scale;
	double V;
	int Nbins, exit_code;
	int *hist;
	double cutoff;
	double *radii, *g;
};

/*
Will check the required_keys have been set.

Inputs:
	* Py_File_Data <PyDictObject*> => The object parsed from the python input.
	* required_keys <*char[]> => The keys that need to be found
	* key_descriptions <*char[]> => The description of each key
Outputs:
	* <int> The number of missing keys
*/
int check_required_keys(PyDictObject* Py_File_Data, char *required_keys[],
					    char *key_descripts[], int num_keys)
{
	PyObject *key, *val;

	// Loop over the required keys and check we have them all
    int keys_missing=0;
    for (int i=0; i<num_keys; i++) {

        key = PyUnicode_FromString(required_keys[i]);

        // Check the key is definitely there.
        if (!PyDict_Contains(Py_File_Data, key)) {
        	printf("\n------------------------------\n");
            printf("\nKey '%s' required in the input dictionary!\n", required_keys[i]);
            printf("\n%s\n\n", key_descripts[i]);
        	printf("\n------------------------------\n");
            keys_missing++;
        }
    }

    return keys_missing;
}

/*
Will create the file_data struct from the input dictionary passed from python.

This function won't do many checks to see if the data looks right, though this API
will be used by the main data analysis code so hopefully other people won't need to
use it.

Inputs:
	* Py_File_Data <PyDictObject*> => The object parsed from the python input.
Outputs:
	* <struct Lammps_Dump> The struct to be used in the RDF C function.
*/
static struct Lammps_Dump set_Lammps_Dump(PyDictObject* Py_File_Data) {
	PyObject *key, *val;
	char *required_keys[] = {"pos", "ABC", "mol_nums", "atom_types"};
	char *key_descripts[] = {
							 "The positions. These should be given as a 2D list of shape (natom, 3).",
							 "The vector that define the simulation box\nThe input should be of the form:\n[    [xlo, xhi, xy],\n     [ylo, yhi, xz],\n     [zlo, zhi, yz] ]",
							 "The molecular number of each atom should be 1D list of length (natom, )",
							 "The type of each atom, should be 1D list of length (natom, )",
							};

	struct Lammps_Dump file_data;
	file_data.exit_code = 0;

    // Loop over the required keys and check we have them all
    int num_keys = (int) sizeof(required_keys) / sizeof(required_keys[0]);
    int keys_missing = check_required_keys(Py_File_Data, required_keys, key_descripts, num_keys);

    // Return an error is required
    if (keys_missing > 0) {
      printf("\n\n%d keys missing from the input dictionary!\n\n", keys_missing);
      file_data.exit_code = 1;
      return file_data;
    }

    /*
    Setting data within the Lammps_Dump struct
    */

    // Get natom
    key = PyUnicode_FromString("pos");
    val = PyDict_GetItem(Py_File_Data, key);
    file_data.natom = PyObject_Length(val);

    // Get timestep
    key = PyUnicode_FromString("timestep");
    if (PyDict_Contains(Py_File_Data, key)) {
      val = PyDict_GetItem(Py_File_Data, key);
      file_data.timestep = (int) PyLong_AsLong(val);
    }

    // Get is_triclinic
    key = PyUnicode_FromString("is_triclinic");
    if (PyDict_Contains(Py_File_Data, key)) {
      val = PyDict_GetItem(Py_File_Data, key);
      file_data.is_triclinic = (int) PyLong_AsLong(val);
    } else { file_data.is_triclinic = 1; }

    // Get xlo, xhi, ylo, yhi etc...
    // ABC = [[xlo, xhi, xy],
    //        [ylo, yhi, xz],
    //        [zlo, zhi, yz] ]
    key = PyUnicode_FromString("ABC");
    val = PyDict_GetItem(Py_File_Data, key);
    int len_ABC = PyObject_Length(val);
    if (len_ABC != 3) {
    	printf("\nThe length of your ABC array isn't correct.\n\n");
    	printf("\nThe input should be of the form:\n");
    	printf("\n[    [xlo, xhi, xy],");
    	printf("\n     [ylo, yhi, xz],");
    	printf("\n     [zlo, zhi, yz] ]\n");
    	file_data.exit_code = 2;
    	return file_data;
    }
    file_data.xlo = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 0), 0));
    file_data.xhi = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 0), 1));
    file_data.xy = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 0), 2));
    file_data.ylo = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 1), 0));
    file_data.yhi = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 1), 1));
    file_data.xz = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 1), 2));
    file_data.zlo = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 2), 0));
    file_data.zhi = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 2), 1));
    file_data.yz = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, 2), 2));

    // Allocate the arrays
    file_data.R.x = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.R.y = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.R.z = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.sR.x = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.sR.y = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.sR.z = (double*) malloc ( file_data.natom * sizeof(double) );
    file_data.mol = (int*) malloc( file_data.natom * sizeof(int) );
    file_data.type = (int*) malloc( file_data.natom * sizeof(int) );

    // Check the lengths of all the arrays
    key = PyUnicode_FromString("atom_types");
    val = PyDict_GetItem(Py_File_Data, key);
    int len_types = PyObject_Length(val);
    key = PyUnicode_FromString("mol_nums");
    val = PyDict_GetItem(Py_File_Data, key);
    int len_mols = PyObject_Length(val);
    if (file_data.natom != len_mols || file_data.natom != len_types) {
      printf("\n\nYour arrays are not the same length!\n\n");
      printf("\n\nPlease make sure they are the same length when passed into the calc_RDF func.");
      printf("\n    len(pos)        = %d\n    len(mol_nums)   = %d\n    len(atom_types) = %d\n\n",
             file_data.natom, len_mols, len_types);
      file_data.exit_code = 3;
      return file_data;
    }
    // Set the positions
    key = PyUnicode_FromString("pos");
    val = PyDict_GetItem(Py_File_Data, key);
    for (int i=0; i<file_data.natom; i++) {
        file_data.R.x[i] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, i), 0));
        file_data.R.y[i] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, i), 1));
        file_data.R.z[i] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(val, i), 2));
    }

    key = PyUnicode_FromString("atom_types");
    val = PyDict_GetItem(Py_File_Data, key);
    for (int i=0; i<file_data.natom; i++) {
        file_data.type[i] = (int) PyLong_AsLong(PyList_GetItem(val, i));
    }
    key = PyUnicode_FromString("mol_nums");
    val = PyDict_GetItem(Py_File_Data, key);
    for (int i=0; i<file_data.natom; i++) {
        file_data.mol[i] = (int) PyLong_AsLong(PyList_GetItem(val, i));
    }

	return file_data;
}

/*
Will set the params struct

This function won't do many checks to see if the data looks right, though this API
will be used by the main data analysis code so hopefully other people won't need to
use it.

Inputs:
	* Py_File_Data <PyDictObject*> => The object parsed from the python input.
Outputs:
	* <struct Params> The struct to be used in the RDF C function.
*/
struct Params set_params(PyDictObject* Py_File_Data) {
	struct Params params;
	PyObject *key, *val;
	char *required_keys[] = {"cutoff", "dr"};
	char *key_descripts[] = {
							 "The cutoff distance used when calculating the RDF (reducing won't speed up the code)",
							 "The spacing between bins in the histogram (increasing won't speed up the code)",
							 };

    // Loop over the required keys and check we have them all
    int num_keys = (int) sizeof(required_keys) / sizeof(required_keys[0]);
    int keys_missing = check_required_keys(Py_File_Data, required_keys, key_descripts, num_keys);

    // Return error if required
    if (keys_missing > 0) {
      printf("\n\n%d keys missing from the input dictionary!\n\n", keys_missing);
      params.exit_code = 1;
      return params;
    }

    // Set the bin spacing
    key = PyUnicode_FromString("dr");
    val = PyDict_GetItem(Py_File_Data, key);
    params.dr = PyFloat_AsDouble(val);

    // Set the cutoff
    key = PyUnicode_FromString("cutoff");
    val = PyDict_GetItem(Py_File_Data, key);
    params.cutoff = PyFloat_AsDouble(val);

    // Get scale
    key = PyUnicode_FromString("scale");
    if (PyDict_Contains(Py_File_Data, key)) {
      val = PyDict_GetItem(Py_File_Data, key);
      params.scale = PyFloat_AsDouble(val);
    } else { params.scale = 1.0; }

	return params;
}


/*
Will set the config file struct

This function won't do many checks to see if the data looks right, though this API
will be used by the main data analysis code so hopefully other people won't need to
use it.

Inputs:
	* Py_File_Data <PyDictObject*> => The object parsed from the python input.
Outputs:
	* <struct Config_File> The struct to be used in the RDF C function.
*/
struct Config_File set_config_struct(PyDictObject* Py_File_Data) {
	struct Config_File conf_file;
	PyObject *key, *val;
	int len_types;
	char *required_keys[] = {"type_list1", "type_list2"};
	char *key_descripts[] = {
							 "The atom types to calculate the RDF for",
							 "The atom types to calculate the RDF for",
							 };

    // Loop over the required keys and check we have them all
    int num_keys = (int) sizeof(required_keys) / sizeof(required_keys[0]);
    int keys_missing = check_required_keys(Py_File_Data, required_keys, key_descripts, num_keys);

    // Return error if required
    if (keys_missing > 0) {
      printf("\n\n%d keys missing from the input dictionary!\n\n", keys_missing);
      conf_file.exit_code = 1;
      return conf_file;
    }

    // Set the first type list
    key = PyUnicode_FromString("type_list1");
    val = PyDict_GetItem(Py_File_Data, key);
    len_types = PyObject_Length(val);
    conf_file.type_1.N = len_types;
    conf_file.type_1.list = (int*) malloc( len_types * sizeof(int) );
    for (int i=0; i<len_types; i++) {
    	conf_file.type_1.list[i] = (int) PyLong_AsLong(PyList_GetItem(val, i));
    }

    // // Set the second type list
    key = PyUnicode_FromString("type_list2");
    val = PyDict_GetItem(Py_File_Data, key);
    len_types = PyObject_Length(val);
    conf_file.type_2.N = len_types;
    conf_file.type_2.list = (int*) malloc( len_types * sizeof(int) );
    for (int i=0; i<len_types; i++) {
    	conf_file.type_2.list[i] = (int) PyLong_AsLong(PyList_GetItem(val, i));
    }

	return conf_file;
}

/*
    Calculate cos(angle) between vecs
*/
double proj(struct Single_Pos a, struct Single_Pos b);
double proj(struct Single_Pos a, struct Single_Pos b)
{
    double dot,norm;
    dot =  (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    norm = (a.x * a.x) + (a.y * a.y) + (a.z * a.z);
    return dot / norm;
}


/*
Will normalise a vector.

Inputs:
    * vec <struct Single_Pos*> => The vector struct to normalise.
Outputs:
    <struct Single_Pos*> normalised vector
*/
struct Single_Pos normalise_vec( struct Single_Pos vec) {
    double norm = sqrt( (vec.x*vec.x) + (vec.y*vec.y) + (vec.z*vec.z) );
    
    struct Single_Pos normed_vec;

    normed_vec.x = vec.x / norm;
    normed_vec.y = vec.y / norm;
    normed_vec.z = vec.z / norm;

    return normed_vec;
}

/*
Will calculate the cell vectors from the dump data.

This will take the cell vector info from the lammps dump file
and calculate the ABC vectors for the cell with it.

Inputs
    * file_data <struct Lammps_Dump*> => The data from the file.
Outputs
    <struct Matrix> The ABC cell vectors.
*/
struct Matrix get_cell_vecs(struct Lammps_Dump *file_data) {
    struct Matrix ABC;

    // Get oblique supercell vectors
    ABC.a.x = file_data->xhi - file_data->xlo;
    ABC.a.y = 0.0; ABC.a.z = 0.0; 

    ABC.b.x = file_data->xy;
    ABC.b.y = file_data->yhi - file_data->ylo;
    ABC.b.z = 0.0;

    ABC.c.x = file_data->xz;
    ABC.c.y = file_data->yz;
    ABC.c.z = file_data->zhi - file_data->zlo;

    return ABC;
}

/*
Will check if a file exists

Inputs:
    * filepath <char[]> => The path to the file that needs checking.
Outputs:
    <int> 1 if the file exists, 0 if not.
*/
int file_exists (char filepath[]) {
    FILE *file_obj = fopen(filepath, "r");

    if (file_obj == NULL) {
        return 0;
    } else {
        fclose(file_obj);
        return 1;
    }
}

/*
    Will read a config file.

    The config file has parameters telling the code what atom types we want to find the rdf from.

    It has 2 sections: type_1 and type_2.
    
    In each section we declare a number on each line. The first number gives the number of types
    we would like to declare. On the next N lines we should give the types that we would like to
    declare.

    The format is something like:
    `
    type_1
    1     # num types
    1     # type 1
    type_2
    2     # 2 types
    1     # type 1
    2     # type 2
    `
    This would mean that we have 1 type (1) in type_1 and 2 types (1, 2) in type_2.

    Inputs:
        * file_path <char[]> => The path to the config file
        * config_file <Config_File*> => The struct holding all the data from the file.
    Outputs:
        <int> exit_code... Non-zero means there's been an error.
*/
struct Config_File read_config ( char file_path[]) {

    FILE *file_obj;
    struct Config_File config_file;
    char buffer[cmax_length], *buff2;
    int i = 0;

    config_file.exit_code = 0;

    // Some error checking and open file.
    if (file_exists(file_path) == 0) {
        printf("\n\nCan't find file: '%s'\n", file_path);
        config_file.exit_code = 1;
        return config_file;
    }
    file_obj = fopen(file_path, "r");

    /*
        Get info for type 1 
    */
    // First find the start of type_1
    int found_type = 0;
    while( fgets(buffer, cmax_length, file_obj) != NULL) {
        if ( strcmp(buffer, "type_1\n") == 0 ) {
            found_type = 1;
            break;
        }
    }
    // Check we found the required section.
    if (found_type == 0) { 
        printf("\n\nCan't find type_1 section\n");
        config_file.exit_code = 2;
        fclose(file_obj);
        return config_file;
    }

    // Get number of types
    buff2 = fgets(buffer, cmax_length, file_obj);
    sscanf(buffer,"%d",&config_file.type_1.N);
    config_file.type_1.list = (int*) malloc(config_file.type_1.N * sizeof(int));


    // Loop over the next N lines to get each type.
    for(i=0; i < config_file.type_1.N; ++i)
    {
        buff2 = fgets(buffer, cmax_length, file_obj);
        sscanf(buffer, "%d", &config_file.type_1.list[i]);
    }
    rewind(file_obj);

    /*
        Get info for type 2
    */
    // Find start of type_2 section
    found_type = 0;
    while( fgets(buffer, cmax_length, file_obj) != NULL) {
        if( strcmp(buffer,"type_2\n") == 0 ) {
            found_type = 1;
            break;
        }
    }
    // Check we found the required section.
    if (found_type == 0) { 
        printf("\n\nCan't find type_2 section\n");
        config_file.exit_code = 2;
        fclose(file_obj);
        return config_file;
    }

    buff2 = fgets(buffer, cmax_length, file_obj);
    sscanf(buffer, "%d", &config_file.type_2.N);
    config_file.type_2.list = (int*) malloc(config_file.type_2.N * sizeof(int));
    for(i = 0;i<=config_file.type_2.N;++i)
    {
        buff2 = fgets(buffer, cmax_length, file_obj);
        sscanf(buffer, "%d", &config_file.type_2.list[i]);
    }
    fclose(file_obj);

    return config_file;
}

/*
Will read a lammps file from a filepath and return a Lammps_Dump object

This will simply iterate over all lines and call buff2 = fgets on each one to parse
them.
The unit cell vectors will also be calculated and stored in the Lammps_Dump
struct.

Nothing is returned, the second (pointer) argument, file_data, is passed as 
a pointer and modified within the function

Inputs:
    * filepath <char[]> => The path pointing towards the lammps file.
    * file_data <struct Lammps_Dump*> => The struct object to hold the file data.
Outputs:
    <Lammps_Dump> The struct holding all the lammps file file_data.
*/
struct Lammps_Dump read_lammps_dump(char filepath[]) {
    char buffer[cmax_length], *buff2;
    double S[4], min_x, max_x, min_y, max_y;
    int int_buffer, i, nx, ny, nz;
    struct Lammps_Dump file_data;
    FILE *Fin;
    file_data.exit_code = 0;

    // Create the file object
    Fin = fopen(filepath, "r");

    //Check the file exists
    if (Fin == NULL) {
        printf("\n\nCan't find file: '%s'\n", filepath);
        file_data.exit_code = 1;
        fclose(Fin);
        return file_data;
    }

    /*
        Loop over all lines and extract info from each one.

        This is saved in the file_data struct (a Lammps_Dump object).
    */
    // Get step num and natom
    buff2 = fgets(buffer, cmax_length, Fin);
    buff2 = fgets(buffer, cmax_length, Fin);
    sscanf(buffer, "%d", &file_data.timestep);
    buff2 = fgets(buffer, cmax_length, Fin);
    buff2 = fgets(buffer, cmax_length, Fin);
    sscanf(buffer, "%d", &file_data.natom);

    // Get cell vec parameters
    buff2 = fgets(buffer, cmax_length, Fin);
    buff2 = fgets(buffer, cmax_length, Fin);
    sscanf(buffer, "%lf\t%lf\t%lf", &file_data.xlo_bound, &file_data.xhi_bound, &file_data.xy);
    buff2 = fgets(buffer, cmax_length, Fin);
    sscanf(buffer, "%lf\t%lf\t%lf", &file_data.ylo_bound, &file_data.yhi_bound, &file_data.xz);
    buff2 = fgets(buffer, cmax_length, Fin);
    sscanf(buffer, "%lf\t%lf\t%lf", &file_data.zlo_bound, &file_data.zhi_bound, &file_data.yz);

    // If there are no xy, xz, yz vecs then our system isn't triclinic
    if (file_data.xy == 0) {
        file_data.is_triclinic = 0;
    }

    // Get minimum and maximum x
    S[0] = 0.0;
    S[1] = file_data.xy;
    S[2] = file_data.xz;
    S[3] = file_data.xy + file_data.xz;

    min_x = S[0];
    max_x = S[0];
    for (int i=1; i<4; i++) {
        if (S[i] < min_x) { min_x = S[i]; }
        if (S[i] > max_x) { max_x = S[i]; }
    }

    if(file_data.yz<0.0) { min_y = file_data.yz; } else { min_y = 0.0; }
    if(file_data.yz>0.0) { max_y = file_data.yz; } else { max_y = 0.0; }

    file_data.xlo = file_data.xlo_bound - min_x;
    file_data.xhi = file_data.xhi_bound - max_x;
    file_data.ylo = file_data.ylo_bound - min_y;
    file_data.yhi = file_data.yhi_bound - max_y;
    file_data.zlo = file_data.zlo_bound;
    file_data.zhi = file_data.zhi_bound;

    // Allocate memory
    // All positions
    file_data.R.x = (double*) malloc( file_data.natom * sizeof(double) );
    file_data.R.y = (double*) malloc( file_data.natom * sizeof(double) );
    file_data.R.z = (double*) malloc( file_data.natom * sizeof(double) );
    // All scaled positions
    file_data.sR.x = (double*) malloc( file_data.natom * sizeof(double) );
    file_data.sR.y = (double*) malloc( file_data.natom * sizeof(double) );
    file_data.sR.z = (double*) malloc( file_data.natom * sizeof(double) );
    // The molecule number and atom type
    file_data.mol = (int*) malloc( file_data.natom * sizeof(int) );
    file_data.type = (int*) malloc( file_data.natom * sizeof(int) );

    // Read the data
    buff2 = fgets(buffer, cmax_length, Fin);
    for(i=0; i<=file_data.natom; ++i) {
        // read
        buff2 = fgets(buffer, cmax_length, Fin);
        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d",
                      &int_buffer, &file_data.mol[i], &file_data.type[i], 
                      &file_data.R.x[i], &file_data.R.y[i], &file_data.R.z[i],
                      &nx, &ny, &nz);
    } 

    if (fgets(buffer, cmax_length, Fin)) {
        printf("\n\n\nCorrupt Snapshot file. There's more to read at the bottom of the file.\n\n\n");
        printf("\nThis may be OK, but beware of odd results!\n\n");
    }
    fclose(Fin);
    return file_data;
}

/*
Will get the numbers of atoms that belong to the types in the config file.

Inputs:
    * file_data <struct Lammps_Dump*> => The data from Lammps dump file.
    * config_file <struct Config_File*> => The params from the config file
Outputs:
    <struct belongs_to_lists> The number of atoms of each type and which ones.
*/
struct belongs_to_lists get_belongs_to (struct Lammps_Dump *file_data,
                                        struct Config_File *config_file) {

    struct belongs_to_lists belongs_to;
    size_t arr_size = file_data->natom * sizeof(int);

    belongs_to.N1 = 0;
    belongs_to.N2 = 0;
    belongs_to.list1 = (int*)  malloc(arr_size);
    belongs_to.list2 = (int*)  malloc(arr_size);

    // Loop over all atoms
    for (int i=0; i<file_data->natom; i++) {
        // Initialise the arrays
        belongs_to.list1[i] = 0;
        belongs_to.list2[i] = 0;

        // Check if the atom is in any type1 asked about
        for (int j=0; j<config_file->type_1.N; j++) {

            if ( file_data->type[i] == config_file->type_1.list[j] ) {
                belongs_to.N1++;
                belongs_to.list1[i] = 1;
                break;
            }
        }
        // Check if the atom is in any type2 asked about
        for (int j=0; j<config_file->type_2.N; j++) {

            if ( file_data->type[i] == config_file->type_2.list[j] ) {
                belongs_to.N2++;
                belongs_to.list2[i] = 1;
                break;
            }
        }
    }

    return belongs_to;
}

/*
A huge (somewhat ugly) function to wrap the coordinates back into the simulation box.

This will calculate the cell vectors and use them to wrap the coordinates.

Inputs:
    * file_data <struct Lammps_Dump*> => The lammps dump file
Outputs:
    <void> No output, will change the coords within the file_data input.
*/
void wrap_coords(struct Lammps_Dump *file_data) {

    double projaa0, projab0, projac0,
           projba0, projbb0, projbc0,
           projca0, projcb0, projcc0;
    double denom, a_norm_REF, b_norm_REF, c_norm_REF;
    double ih_11, ih_12, ih_13, ih_22, ih_23, ih_33;

    struct Matrix ABC_hat;
    struct Matrix ABC_hat_2;
    struct Matrix Identity;
    struct Matrix ABC = get_cell_vecs(file_data);

    // Get unit vecs
    ABC_hat.a = normalise_vec(ABC.a);
    ABC_hat.b = normalise_vec(ABC.b);
    ABC_hat.c = normalise_vec(ABC.c);

    // Get squared unit vecs for later
    ABC_hat_2.a.x = ABC_hat.a.x * ABC_hat.a.x;
    ABC_hat_2.a.y = ABC_hat.a.y * ABC_hat.a.y;
    ABC_hat_2.a.z = ABC_hat.a.z * ABC_hat.a.z;
    ABC_hat_2.b.x = ABC_hat.b.x * ABC_hat.b.x;
    ABC_hat_2.b.y = ABC_hat.b.y * ABC_hat.b.y;
    ABC_hat_2.b.z = ABC_hat.b.z * ABC_hat.b.z;
    ABC_hat_2.c.x = ABC_hat.c.x * ABC_hat.c.x;
    ABC_hat_2.c.y = ABC_hat.c.y * ABC_hat.c.y;
    ABC_hat_2.c.z = ABC_hat.c.z * ABC_hat.c.z;

    // The identity matrix
    Identity.a.x = 1.0, Identity.a.y = 0.0, Identity.a.z = 0.0;
    Identity.b.x = 0.0, Identity.b.y = 1.0, Identity.b.z = 0.0;
    Identity.c.x = 0.0, Identity.c.y = 0.0, Identity.c.z = 1.0;

    // calculate projections onto the cartesian basis vectors.
    projaa0 = proj(ABC_hat.a, Identity.a);
    projab0 = proj(ABC_hat.a, Identity.b);
    projac0 = proj(ABC_hat.a, Identity.c);
    projba0 = proj(ABC_hat.b, Identity.a);
    projbb0 = proj(ABC_hat.b, Identity.b);
    projbc0 = proj(ABC_hat.b, Identity.c);
    projca0 = proj(ABC_hat.c, Identity.a);
    projcb0 = proj(ABC_hat.c, Identity.b);
    projcc0 = proj(ABC_hat.c, Identity.c);

    // projection denominator (determinant)
    denom = - (projac0 * projbb0 * projca0) \
            + (projab0 * projbc0 * projca0) \
            + (projac0 * projba0 * projcb0) \
            - (projaa0 * projbc0 * projcb0) \
            - (projab0 * projba0 * projcc0) \
            + (projaa0 * projbb0 * projcc0) ;

    // oblique supercell vector lengths
    a_norm_REF = sqrt( (ABC.a.x * ABC.a.x) + (ABC.a.y * ABC.a.y) + (ABC.a.z * ABC.a.z) );
    b_norm_REF = sqrt( (ABC.b.x * ABC.b.x) + (ABC.b.y * ABC.b.y) + (ABC.b.z * ABC.b.z) );
    c_norm_REF = sqrt( (ABC.c.x * ABC.c.x) + (ABC.c.y * ABC.c.y) + (ABC.c.z * ABC.c.z) );

    // inverse box matrix h
    ih_11 = 1.0 / ABC.a.x;
    ih_12 = -ABC.b.x / ( ABC.a.x * ABC.b.y);
    ih_13 = ((ABC.b.x * ABC.c.y) - (ABC.b.y * ABC.c.x)) / (ABC.a.x * ABC.b.y * ABC.c.z);
    ih_22 = 1.0 / ABC.b.y;
    ih_23 = -ABC.c.y / (ABC.b.y * ABC.c.z);
    ih_33 = 1.0 / ABC.c.z;

    // Loop over all atoms.
    for(int j = 0;j<file_data->natom; ++j)
    {
        // projections
        double nom1, nom2, nom3;
        nom1 =  -1.0 * (  (projbc0 * projcb0 * file_data->R.x[j])  \
                        - (projbb0 * projcc0 * file_data->R.x[j])  \
                        - (projbc0 * projca0 * file_data->R.y[j])  \
                        + (projba0 * projcc0 * file_data->R.y[j])  \
                        + (projbb0 * projca0 * file_data->R.z[j])  \
                        - (projba0 * projcb0 * file_data->R.z[j])  );

        nom2 = (  (projac0 * projcb0 * file_data->R.x[j]) \
                - (projab0 * projcc0 * file_data->R.x[j]) \
                - (projac0 * projca0 * file_data->R.y[j]) \
                + (projaa0 * projcc0 * file_data->R.y[j]) \
                + (projab0 * projca0 * file_data->R.z[j]) \
                - (projaa0 * projcb0 * file_data->R.z[j]) );

        nom3 =  - 1.0 * (  (projac0 * projbb0 * file_data->R.x[j]) \
                         - (projab0 * projbc0 * file_data->R.x[j]) \
                         - (projac0 * projba0 * file_data->R.y[j]) \
                         + (projaa0 * projbc0 * file_data->R.y[j]) \
                         + (projab0 * projba0 * file_data->R.z[j]) \
                         - (projaa0 * projbb0 * file_data->R.z[j]) );

        // coordinates with respect to oblique system
        double a, b, c;
        a = nom1 / denom;
        b = nom2 / denom;
        c = nom3 / denom;

        // Squared values
        double a2, b2, c2;
        a2 = a*a; b2 = b*b; c2 = c*c;

        //  + / -  direction
        int a_sign, b_sign, c_sign;
        if( (a * ABC_hat.a.x * ABC.a.x) < 0.0 ) { a_sign =  -1.0; } else { a_sign =  +1.0;}
        if( (b * ABC_hat.b.x * ABC.b.x) + (b * ABC_hat.b.y * ABC.b.y) < 0.0 ) { b_sign =  -1.0; } else { b_sign =  +1.0; }
        if( (c * ABC_hat.c.x * ABC.c.x) + (c * ABC_hat.c.y * ABC.c.y) + (c * ABC_hat.c.z * ABC.c.z) < 0.0 ) { c_sign =  -1.0; } else { c_sign =  +1.0; }

        // image calculation
        double a_norm, b_norm, c_norm;
        a_norm = sqrt( (a2 * ABC_hat_2.a.x) + (a2 * ABC_hat_2.a.y) + (a2 * ABC_hat_2.a.z) );
        b_norm = sqrt( (b2 * ABC_hat_2.b.x) + (b2 * ABC_hat_2.b.y) + (b2 * ABC_hat_2.b.z) );
        c_norm = sqrt( (c2 * ABC_hat_2.c.x) + (c2 * ABC_hat_2.c.y) + (c2 * ABC_hat_2.c.z) );

        file_data->nx = (int) floor( a_norm / a_norm_REF );
        if (a_sign < 0.0) file_data->nx = -file_data->nx - 1;
        file_data->ny = (int) floor( b_norm / b_norm_REF );
        if (b_sign < 0.0) file_data->ny = -file_data->ny - 1;
        int nz = (int) floor( c_norm / c_norm_REF );
        if( c_sign < 0.0) nz =  -nz - 1;

        // periodic wrapping
        file_data->R.x[j] = file_data->R.x[j] - (file_data->nx * ABC.a.x) - (file_data->nz * ABC.b.x) - (nz * ABC.c.x);
        file_data->R.y[j] = file_data->R.y[j] - (file_data->nz * ABC.b.y) - (nz * ABC.c.y);
        file_data->R.z[j] = file_data->R.z[j] - (nz * ABC.c.z);

        // scaled coordinates
        file_data->sR.x[j] = (ih_11 * file_data->R.x[j]) + (ih_12 * file_data->R.y[j]) + (ih_13 * file_data->R.z[j]);
        file_data->sR.y[j] = (ih_22 * file_data->R.y[j]) + (ih_23 * file_data->R.z[j]);
        file_data->sR.z[j] = ih_33 * file_data->R.z[j];
    }
}

/*
    Will get the histogram count.

    This will wrap the coordinates first, using the scaled vectors, and
    then use these scaled vectors to find the distances to build the
    histogram.
    
    Inputs:
        * file_data <Lammps_Dump> => The file_data from the snapshot file
        * bin <int> => The number of bins
    Outputs:
        <int*> The bin counts
*/
void get_histogram(struct Lammps_Dump file_data, int bins, int *counts, double dr, struct belongs_to_lists belongs_to) {
    int diff_mol, correct_atom_type, bh, i, j;
    double dist, rx, ry, rz;
    double svectorx, svectory, svectorz;
    struct Matrix ABC = get_cell_vecs(&file_data);

    // Init the histogram at the beginning of each step
    for (i = 0; i<bins; ++i){ counts[i] = 0; }

    // Loop over all pairs of atoms
    for(j = 0; j<file_data.natom-1; ++j) {
        printf("\rAtom %d/%d                   \r", j, file_data.natom - 1);

        for(i=j+1; i<file_data.natom; ++i) {

            // If the pair of atoms are the correct type and on different mols then count them
            diff_mol = file_data.mol[i] != file_data.mol[j];
            correct_atom_type = (belongs_to.list1[i] && belongs_to.list2[j]) || (belongs_to.list1[j] && belongs_to.list2[i]);
            if(diff_mol && correct_atom_type) {

                // svectors are for the coordinate wrapping
                svectorx = file_data.sR.x[j] - file_data.sR.x[i];
                svectory = file_data.sR.y[j] - file_data.sR.y[i];
                svectorz = file_data.sR.z[j] - file_data.sR.z[i];

                svectorx = svectorx - rint(svectorx);
                svectory = svectory - rint(svectory);
                svectorz = svectorz - rint(svectorz);

                // The actual displacement vectors
                rx = (ABC.a.x * svectorx) + (file_data.xy * svectory) + (file_data.xz * svectorz);
                ry = (ABC.b.y * svectory) + (file_data.yz * svectorz);
                rz = (ABC.c.z * svectorz);

                dist = sqrt( (rx*rx) + (ry*ry) + (rz*rz) );

                // Add to the count
                bh = floor(dist / dr);
                if (bh < bins) { counts[bh] = counts[bh] + 2; }
            }
        }
    }
}

/*
Will calculate the RDF from the inputted data.

Inputs:
    * file_data <struct Lammps_Dump> => The dump data from the lammps snapshot file.
    * config_file <struct Config_File> => The config file data.
    * params <struct Params*> => The parameters used to calculate the RDF.
Outputs:
    <double *> The RDF.
*/
int calc_RDF(struct Lammps_Dump file_data, struct Config_File config_file, struct Params *params) {

    // number of bins
    params->Nbins = floor( (params->cutoff / params->dr) + 0.5);
    params->g = (double*)  malloc( params->Nbins * sizeof(double) );
    params->radii = (double*)  malloc( params->Nbins*sizeof(double) );
    params->hist = (int*) malloc( params->Nbins * sizeof(int) );
    

    for (int i=0; i<params->Nbins; ++i) {
        params->radii[i] = (i + 0.5) * params->dr;
        params->g[i] = 0.0;
    }

    // Shift so xmin, ymin, zmin are at the orign.
    for(int j = 0;j<file_data.natom;++j)
    {
        file_data.R.x[j] = file_data.R.x[j] - file_data.xlo;
        file_data.R.y[j] = file_data.R.y[j] - file_data.ylo;
        file_data.R.z[j] = file_data.R.z[j] - file_data.zlo;
    }

    // Get volume
    struct Matrix ABC = get_cell_vecs(&file_data);
    params->V = ABC.a.x * ABC.b.y * ABC.c.z;

    // Wrap coordinates
    wrap_coords(&file_data);

    // Get which atoms belong to which type and how many
    struct belongs_to_lists belongs_to = get_belongs_to( &file_data, &config_file );

    // Get histogram
    get_histogram(file_data, params->Nbins, params->hist, params->dr, belongs_to);

    // rdf calculation
    double norm = params->V / (4.0 * pi * belongs_to.N1 * belongs_to.N2 * params->scale * params->dr);
    for(int i=0; i<params->Nbins; i++) {
        params->g[i] = params->g[i] + (params->hist[i] * norm) / (params->radii[i] * params->radii[i]);
    }

    // Free the belongs to data
    free(belongs_to.list1);
    free(belongs_to.list2);
    // Free the Lammps Data
    free(file_data.R.x);
    free(file_data.R.y);
    free(file_data.R.z);
    free(file_data.mol);
    free(file_data.type);
    free(file_data.sR.x);
    free(file_data.sR.y);
    free(file_data.sR.z);

    return 0;
}
