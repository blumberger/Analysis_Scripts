#include <Python.h>
#include <stdio.h>
#include "RDF_utils.h"


// Will accept a dict as an argument and print the key x
static PyObject *
wrap_coords_module_wrap_coords(PyObject *self, PyObject *args)
{
    PyDictObject* Py_File_Data;

    // Parse the arguments
    if (!PyArg_ParseTuple(args, "O", &Py_File_Data)) {
        return NULL;
    }


    // Set the struct that contains all the position data.
    struct Lammps_Dump file_data = set_Lammps_Dump(Py_File_Data);
    if (file_data.exit_code != 0) {
      printf("\n\nError setting file_data: Exitting Gracefully\n\n");
      return NULL;   }


    // PyObject  *list;
    // PyObject *rdf = PyList_New(params.Nbins);
    // for (int i=0; i<params.Nbins; i++) {
    //     list = Py_BuildValue("f", params.g[i]);
    //     PyList_SetItem(rdf, i, list);
    // }
    // PyObject *radii = PyList_New(params.Nbins);
    // for (int i=0; i<params.Nbins; i++) {
    //     list = Py_BuildValue("f", params.radii[i]);
    //     PyList_SetItem(radii, i, list);
    // }
    // PyObject *return_list = PyList_New(2);
    // PyList_SetItem(return_list, 0, rdf);
    // PyList_SetItem(return_list, 1, radii);

    // return return_list;
    Py_RETURN_NONE;
}


//Method definition object for this extension, these argumens mean:
//ml_name: The name of the method
//ml_meth: Function pointer to the method implementation
//ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
//ml_doc:  Contents of this method's docstring
static PyMethodDef RDF_module_methods[] = { 
    {
         "wrap_coords",
         wrap_coords_module_wrap_coords,
         METH_VARARGS,
         "Wrap triclinic coords."
    },
    {NULL, NULL, 0, NULL}
};

//Module definition
//The arguments of this structure tell Python what to call your extension,
//what it's methods are and where to look for it's method definitions
static struct PyModuleDef RDF_module_definition = { 
    PyModuleDef_HEAD_INIT,
    "RDF_module",
    "A Python module to unwrap coordinates to create whole molecules again"
    -1, 
    RDF_module_methods
};

//Module initialization
//Python calls this function when importing your extension. It is important
//that this function is named PyInit_[[your_module_name]] exactly, and matches
//the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_RDF(void)
{
    Py_Initialize();

    return PyModule_Create(&RDF_module_definition);
}
