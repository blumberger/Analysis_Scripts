# How to contribute to the code

In order to make the code long lasting and usable/editable for everyone we need to stick to some ground rules. These are there to ensure the code can be used easily by anyone and that the code doesn't become a sprawling mess.

## The rules

 - **Use an input file**. In order to standardise input for all commands and input file will be used.
 - **Document as you go**. This doesn't just mean comment your code but also: use docstrings; use error checking to print (useful) error messages to let the user know what is wrong with their input file; save well explained example input files for any new features you implement.
 - **Use clear variable names**. Using clear variable names helps a lot when trying to understand what the script is doing. They can often replace comments (which saves time and effort for the developer)! See [PEP8](https://realpython.com/python-pep8/#naming-styles) for more info. P.S. this extends to function names, filenames, class names etc...
 - **Write with modularity in mind**. Divide up code into sensibly named files. Use objects and inherit from base classes to make common tasks standardised. E.g. calculations should be written in classes derived from the Calc_Type, this has standard names for the data (self.data) etc... that can be accessed later without having to use if statements for special cases.
 - **If the code can be reused, write a function**. If someone wants to write an xyz file and the function that writes them is incorporated into another function they have to do some tedious editting. Remember to always give that function a docstring though! They really help.
 - **Follow PEP8**. If you are writing in python, try and follow PEP8 guidelines as closely as possible. Most IDEs have PEP8 style checking, turn it on and try to follow it. Of course, if you think for a special case PEP8 obfuscates the meaning don't use it, though the default should be PEP8. Also don't use start imports (i.e. from ... import \*). Finally remember whitespace can make things so much nicer to read.
 - **If your writing in python: Use python3**. Any development on Python2 has stopped and the langauge has become deprecated. It will still work as it should do but won't get better. Python3 also has lots of nice new features.




## Where to put new code?
Put it in the src directory in an appropriately named folder. All of the code will be called from the main.py file so make sure it can be imported from this file.

## How to document your new code?
As you add functionality you should test it in your own input file. E.g. if you add some code to load and write .dat files then write an input file that does that to run t. After you are happy your code and the input file commands work comment the input file well and copy it into an appropriate directory in the examples directory.

## What should go in the repo?
Any code you think may be helpful, can be ran with an input file and is documented.
