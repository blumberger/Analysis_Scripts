import os
import subprocess
import re

from src.parsing import general_parsing  as gen_parse

def test_out_files():
    if os.path.isfile("ammended.xyz"):
        os.remove
        return True

true = lambda std_out: True

def check_stdout(std_out):
    """
    Will check the std_out to see if everything has been printed as expected.
    """
    ltxt = std_out.split("\n")
    errs = []
    for line in ltxt:
        should_be = re.findall("\(.*should *be.*\)", line, flags=re.IGNORECASE)
        if len(should_be) != 1:
            continue

        should_be = should_be[0]
        out = line.replace(should_be, "").strip()

        should_be = re.findall("d *be.*\)", should_be)[0][:-1]
        should_be = should_be[should_be.find("be")+2:].strip()
        should_be = gen_parse.rm_quotation_marks(should_be)

        if out != should_be:
            errs.append(f"Line might be incorrect: {line}." +
                        "\n\nOutput is" + f"{out} though it seems like it "
                        + "should be {should_be}.")

    return errs


# Set up the test folders
examples_folder = './examples'


# Ammend this dictionary to add new tests:
#    Keys  = foldername within examples dir
#    Value = dict with keys:
#        test_func = a function to call to test (use true func if you can't think of one!)
#        inp_file  = the filename of the input file to run
#        out_file  = the filename of the output file that is created (use False if there isn't one)
all_example_tests = {
                        'ammending_xyz_files': {
                                        'test_func': check_stdout,
                                        'inp_file': 'ammend_xyz_file.inp',
                                        'out_file': 'ammended.xyz',
                                               },

                        'calc_pvecs': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_pvecs.inp',
                                        'out_file' : 'pvecs.xyz',
                                      },

                        'Nearest_Neighbour':   {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_NN.inp',
                                        'out_file' : 'NN.json',
                                               },
                        'prettify_CP2K_input': {
                                        'test_func':  check_stdout,
                                        'inp_file' : 'prettify_CP2K_input.inp',
                                        'out_file' : 'pretty.inp',
                                               },

                        'arithmetic': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'arithmetic.inp',
                                        'out_file' : False,
                                      },
                        'metadata':   {
                                        'test_func': check_stdout,
                                        'inp_file' : 'metadata_setting.inp',
                                        'out_file' : False,
                                      },
                        'calc_density': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_densities.inp',
                                        'out_file' : 'density.csv',
                                        },
                        'unwrap_snapshot_data': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'unwrap_snapshot_data.inp',
                                        'out_file' : ('unwrapped.xyz', "wrapped.xyz"),
                                                },
                      'calc_angular_distribution': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_angular_dist.inp',
                                        'out_file' : 'ang_dist.json',
                                                   },
                        'for_loops': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'for_loops.inp',
                                        'out_file' : False,
                                     },
                        'calc_RDF': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_RDF.inp',
                                        'out_file' : "rdf.csv",
                                    },

                        'scripts': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'run_script.inp',
                                        'out_file' : False,
                                   },
                        'create_psf_file_from_xyz':  {
                                        'test_func': check_stdout,
                                        'inp_file': 'create_psf_from_xyz.inp',
                                        'out_file': ("pentacene_charged.psf", "pentacene_neutral.psf", "bond_struct.png"),
                                                     },
                        'inline_python_script': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'run.inp',
                                        'out_file' : False,
                                                },
                        'calc_coupling_connections': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'calc_couplings.inp',
                                        'out_file' : [f"Hab_connect_layer_{i}.png" for i in range(0, 8)],
                                                     },     
                        'calc_layers': {
                                        'test_func': check_stdout,
                                        'inp_file' : 'get_layers.inp',
                                        'out_file' : False,
                                       },


                      # Keep this one last
                      'ipython_shell':      {
                                             'test_func': check_stdout,
                                             'inp_file' : 'IPython_shell.inp',
                                             'out_file' : False,
                                             'in_shell' : True,
                                            },
                    }

all_example_tests = {f'{examples_folder}/{i}': all_example_tests[i] for i in all_example_tests}
all_errors = {}

# Loop over the example folders and run them all
for test_fold in all_example_tests:
    print("#################################################")
    print(f"  Testing folder '{test_fold.split('/')[-1]}'")

    # Set some initial variables
    test_dict = all_example_tests[test_fold]
    test_func = test_dict['test_func']
    inp_file = f"{test_fold}/{test_dict['inp_file']}"
    out_file = test_dict['out_file']
    interactive = test_dict.setdefault('in_shell', False)

    # Check the inp file is there
    if not os.path.isfile(inp_file):
        all_errors.setdefault(test_fold, []).append(f"Couldn't find an input file in the folder {test_fold}. Looked for {inp_file}")
        continue

    # Run the input file
    print(" ")
    cmd = ["bash", "run.sh", "-i", f"{inp_file}"]
    print(f"  Calling command '{' '.join(cmd)}'")
    if interactive:
        subprocess.call(cmd)
    else:
        pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        std_out, std_err = pipes.communicate()
        std_err = std_err.decode("utf-8")
        if std_err:
            print("\n\nERROR OCCURED\n\n")
            all_errors[test_fold] = [std_err]

        std_out = std_out.decode("utf-8")

        # Call the test function
        errs = test_func(std_out)
        if errs:
            for err in errs:
                all_errors.setdefault(test_fold, []).append(err)


    # Check the file has been created
    if out_file:
        if type(out_file) == str:
            out_file = (out_file,)
        for file_ in out_file:
            if not os.path.isfile(file_):
                all_errors.setdefault(test_fold, []).append(f"Couldn't find the output file {file_}")
                continue
            else:
                os.remove(file_)
                print(" ")
                print(f"  Removed {file_}")


    print("\n\n#################################################\n\n\n\n\n")


if all_errors:
    print("\n\n\n---\nError Report:\n")
    for test_fold in all_errors:
        print(f"Error {test_fold}:")
        for i, err in enumerate(all_errors[test_fold]):
            print(f"Error {i}: {err}")
else:
    print("All Successful!")
