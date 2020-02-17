# How to contribute to the code

In order to make the code long lasting and usable/editable for everyone we need to stick to some ground rules. These are there to ensure the code can be used easily by anyone and that the code doesn't become a sprawling mess.

# The rules

 - **Use an input file**. In order to standardise input for all commands and input file will be used.
 - **Document as you go**. This doesn't just mean comment your code but also: use docstrings; use error checking to print (useful) error messages to let the user know what is wrong with their input file; save well explained example input files for any new features you implement.
 - **Use clear variable names**. Using clear variable names helps a lot when trying to understand what the script is doing. They can often replace comments (which saves time and effort for the developer)! See [PEP8](https://realpython.com/python-pep8/#naming-styles) for more info. P.S. this extends to function names, filenames, class names etc...
 - **Write with modularity in mind**. Divide up code into sensibly named files. Use objects and inherit from base classes to make common tasks standardised. E.g. calculations should be written in classes derived from the Calc_Type, this has standard names for the data (self.data) etc... that can be accessed later without having to use if statements for special cases.
 - **If the code can be reused, write a function**. If someone wants to write an xyz file and the function that writes them is incorporated into another function they have to do some tedious editting. Remember to always give that function a docstring though! They really help.
 - **Follow PEP8**. If you are writing in python, try and follow PEP8 guidelines as closely as possible. Most IDEs have PEP8 style checking, turn it on and try to follow it. Of course, if you think for a special case PEP8 obfuscates the meaning don't use it, though the default should be PEP8. Also don't use start imports (i.e. from ... import \*). Finally remember whitespace can make things so much nicer to read.
 - **If your writing in python: Use python3**. Any development on Python2 has stopped and the langauge has become deprecated. It will still work as it should do but won't get better. Python3 also has lots of nice new features.


# How the code works
There should be only 2 points of contact that the user has with the code. First they should run the run.sh file. This should point to the second -the input file.

When the run.sh file is called it carries out some initialisation. This involves creating a virtual environment and installing the required python libraries (if it is the first time the code is ran) and then it parses the flags. The available flags are -i and -t. -i points towards the input file. -t signals that we want to run the tests.

If the -t is specified the bash file will call the test files from the folder src/tests.

If the -i is used, then the bash script will then call the main.py script. This main.py serves only to glue everything together. It parses the arguments served to it, imports the input file parser and parses the input file. The code that carries out the input file parsing can be found in src/input_file/input_file.py.

The input file parser will loop over all lines in the input file twice. The first loop is a quick loop to check for errors such as bad syntax or undefined variables. The error messages should be as clear as possible here. To print error messages use the function `self.__print__error()` in the input file parser.

The second loop will parse the input file and carry out any commands. 


# Testing the code after changes
In order to check everything is working before committing the code to the master branch run the testing suite and check the output for anything wrong. Do this by running the command:
```
./run.sh -t
```
This will call the script src/tests/call_all_examples.py which will loop over the example input files specified in the dictionary named: `all_example_tests`.

To add your input file that runs the changes you've made to, append an entry to this dictionary. The key should be the folder name of the input file in the examples directory and the value should be another dictionary with keys `test_func` (a function to run to test the input), `inp_file` (the filename of the inp file) and `out_file` (the name of any files that are created).

# Q & A

## Where to put new code?
Put it in the src directory in an appropriately named folder. All of the code will be called from the main.py file so make sure it can be imported from this file.

## How to document your new code?
As you add functionality you should test it in your own input file. E.g. if you add some code to load an xyz file, adjust it's values, and write it again then write this in an input file and put it in the examples directory. You should also remember to update the test code to include this input file so its functionality can be tested automatically.

## What should go in the repo?
Any code you think may be helpful, can be ran with an input file and is documented.
