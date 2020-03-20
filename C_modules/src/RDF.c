#include <Python.h>
#include <stdio.h>
#include "RDF_utils.h"

//Actual module method definition - this is the code that will be called by
//RDF_module.print_hello_world
static PyObject * 
RDF_module_print_hello_world(PyObject *self, PyObject *args)
{
    printf("\nHello World!\n");
    Py_RETURN_NONE;
}


// Will accept a dict as an argument and print the key x
static PyObject *
RDF_module_calc_RDF(PyObject *self, PyObject *args)
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

    // Set the params struct
    struct Params params = set_params(Py_File_Data);
    if (params.exit_code != 0) {
      printf("\n\nError setting params: Exitting Gracefully\n\n");
      return NULL;    }

    // Set the config file struct
    struct Config_File conf_file = set_config_struct(Py_File_Data);
    if (conf_file.exit_code != 0) {
      printf("\n\nError setting Config File: Exitting Gracefully\n\n");
      return NULL;    }
   
    int exit_code = calc_RDF(file_data, conf_file, &params);
    if (exit_code) {
      printf("\n\nError calculating RDF: Exitting Gracefully\n\n");
      return NULL;
    }


    PyObject  *list;
    PyObject *rdf = PyList_New(params.Nbins);
    for (int i=0; i<params.Nbins; i++) {
        list = Py_BuildValue("f", params.g[i]);
        PyList_SetItem(rdf, i, list);
    }
    PyObject *radii = PyList_New(params.Nbins);
    for (int i=0; i<params.Nbins; i++) {
        list = Py_BuildValue("f", params.radii[i]);
        PyList_SetItem(radii, i, list);
    }
    PyObject *return_list = PyList_New(2);
    PyList_SetItem(return_list, 0, rdf);
    PyList_SetItem(return_list, 1, radii);

    return return_list;
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
        "print_hello_world",
         RDF_module_print_hello_world,
         METH_NOARGS,
         "Print 'hello world' from a method defined in a C extension.",
    },  
    {
         "calc_RDF",
         RDF_module_calc_RDF,
         METH_VARARGS,
         "Calculate the RDF"
    },
    {NULL, NULL, 0, NULL}
};

//Module definition
//The arguments of this structure tell Python what to call your extension,
//what it's methods are and where to look for it's method definitions
static struct PyModuleDef RDF_module_definition = { 
    PyModuleDef_HEAD_INIT,
    "RDF_module",
    "A Python module to calculate the Radial Distribution Function (RDF) of a molecular system",
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
