import os
import subprocess

def test_out_files():
    if os.path.isfile("ammended.xyz"):
        os.remove
        return True

true = lambda : True

# Set up the test folders
examples_folder = './examples'


# Ammend this dictionary to add new tests:
#    Keys  = foldername within examples dir
#    Value = dict with keys:
#        test_func = a function to call to test (use true func if you can't think of one!)
#        inp_file  = the filename of the input file to run
#        out_file  = the filename of the output file that is created (use False if there isn't one)
all_example_tests = {
                     'ammending_xyz_files': {'test_func': true,
                                             'inp_file': 'ammend_xyz_file.inp',
                                             'out_file': 'ammended.xyz',
                                            },

                     'calc_pvecs': {
                                    'test_func': true,
                                    'inp_file' : 'calc_pvecs.inp',
                                    'out_file' : 'pvecs.xyz',
                                   },

                     'Nearest_Neighbour':   {
                                              'test_func': true,
                                              'inp_file' : 'calc_NN.inp',
                                              'out_file' : 'test.json',
                                            },
                     'prettify_CP2K_input': {
                                             'test_func': true,
                                             'inp_file' : 'prettify_CP2K_input.inp',
                                             'out_file' : 'pretty.inp',
                                             },

                      'arithmetic': {
                                        'test_func': true,
                                        'inp_file' : 'arithmetic.inp',
                                        'out_file' : False,
                                    },
                      'metadata':   {
                                        'test_func': true,
                                        'inp_file' : 'metadata_setting.inp',
                                        'out_file' : False,
                                    },
                      'calc_density':       {
                                             'test_func': true,
                                             'inp_file' : 'calc_densities.inp',
                                             'out_file' : 'density_0.csv',
                                            },
                    }

all_example_tests = {f'{examples_folder}/{i}': all_example_tests[i] for i in all_example_tests}
all_errors = []

# Loop over the example folders and run them all
for test_fold in all_example_tests:
    print("#################################################")
    print(f"  Testing folder '{test_fold.split('/')[-1]}'")

    # Set some initial variables
    test_dict = all_example_tests[test_fold]
    test_func = test_dict['test_func']
    inp_file = f"{test_fold}/{test_dict['inp_file']}"
    out_file = test_dict['out_file']

    # Check the inp file is there
    if not os.path.isfile(inp_file):
        all_errors.append(f"Couldn't find an input file in the folder {test_fold}. Looked for {inp_file}")
        continue 

    # Run the input file
    print(" ")
    print(f"  Calling command 'bash run.sh -i {inp_file}")
    subprocess.call(["bash", "run.sh", "-i", f"{inp_file}"])

    # Call the test function
    test_func()

    # Check the file has been created
    if out_file:
        if not os.path.isfile(out_file):
            all_errors.append(f"Couldn't find the output file {out_file}")
            continue 
        else:
            os.remove(out_file)
            print(" ")
            print(f"  Removed {out_file}")
        
        
    print("\n\n#################################################\n\n\n\n\n")


if all_errors:
    print("\n\n\n---\nError Report:\n")
    for i, err in enumerate(all_errors):
        print(f"Error {i}: {err}")
else:
    print("All Successful!")
