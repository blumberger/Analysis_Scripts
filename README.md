BUG REPORT:
   WARNING: Found another bug in the reading of variables -if there is a / in the variable it doesn't read as a string. Things work when in quotation marks though.
   WARNING: New features not fully tested in this commit

# MD analysis scripts

## Mission Statement
The purpose of this repo is to centralise all the scripts that we use as a group so each time someone runs an MD simulation/SH simulation they don't have to re-write what many people have written before them.

One of the main reasons the analysis scripts are re-written is (I think) because it is often easier (or less frustrating) to write the script yourself rather than learn how to use the script that's already available. Therefore, for this repo I think we need some rules about submitting that make the code both usable for people in the future and readable for anyone wanting to edit and extend the scripts.

To make the code easy to use for potential users a few things are vital (these are expanded on in the GUIDE_FOR_DEVELOPERS.md file):
  * A consistent way to interact with the code (input file).
  * Documentation.
  * Useful error messages.


## How to run the code.

  1) Go to the examples and choose a folder you like the sound of.
  2) Copy the contents of the folder to the root director (where this README file is).
  3) Run the file run.sh with an -i flag giving the input file, e.g. `./run.sh -i <input_file>`

  The first time you run the file it may take a while as it will try to install dependencies.

  To run your own analyses create your own input file and run that.

## The input file

Every interaction you have with the code (unless you are developing it) will be through the input file. You can call it anything you want with any extension as long as you tell the code where it can be found via `./run.sh -i <input_file>`

The input file is read and ran line by line (much like normal code). A new line contains a new command, currently there is no facility to merged long lines together.

For example an input file that looked like:

```
x = 2
echo "x = $x"
x = x+1
echo "x = $x"
```

Would output:
```
x = 2
x = 3
```

You should try this for yourself to convince yourself.

### Syntax:
#### Variable Declarations
Variables are declared using the follow syntax:

```
<var name> = <value>
```

Values can be strings, floats or ints. Currently lists are not supported.

Variables can be overwritten by re-declaring them.

#### Variable modification
##### Metadata
Variables can have metadata associated with them and this can be assigned with the following syntax (much like python dictionaries):

```
<var name>['<metadata_name>'] = <value>
```

This is used for some calculations that need specific information that isn't read in via a file.
##### Arithmetic
If the variables declared or loaded from files can be manipulated mathematically then you can do so using normal mathematical expressions e.g:
```
x = 1+2
echo $x
x = x *2
echo $x
x = 3**2
echo $x
x = x/2
```

#### Loading data
Data can be loaded from a file via the following syntax:
```
load <filepath> <file_type> as <var_name>
```
Here filepath points towards the file you would like to load, file_type is the type of file you would like to load (see example for more info) and var_name is the name associated to the data.
