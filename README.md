**BUG REPORT:**

   WARNING: Pvecs calculator not fully tested

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

Examples of input file are given in the examples folder. To run any one of them run the command:
```
./run.sh -i <example input filepath>
```

### A quick example:
An input file that did some basic arithmetic might look like the one below:
```
x = 2
echo "x = $x"
x = x+1
echo "x = $x"
```

This would output:
```
x = 2
x = 3
```
You can try this for yourself by creating an input file with the above text and run it by using the run.sh bash script.

### Syntax:
#### Variable Declarations
Variables are declared using the follow syntax:

```
<variable name> = <value>
```

Values can be strings, floats, ints or lists.

Strings are set either with or without double quotation marks: "str"

floats and ints are set without quotation marks

lists are set with commas i.e. bob, 2, 9 (currently nested lists are not implemented)

Variables can be overwritten by re-declaring them.

Examples can be found in most of the example input files in the examples folder.

#### Variable modification
##### Metadata
Variables can have metadata associated with them and this can be assigned with the following syntax (much like python dictionaries):

```
<variable name>['<metadata_name>'] = <value>
```

This is used for some calculations that need specific information that isn't read in via a file.

An example of metadata being used can be found in the calc_pvecs example input file.

##### Arithmetic
If the variables declared or loaded from files can be manipulated mathematically then you can do so using normal mathematical expressions e.g:
```
x = 1+2
echo $x                  (Should be 3)

x = x * 2
echo $x                  (Should be 6)

x = 3**2                     
echo $x                  (Should be 9)

x = (3 * (7-6) + $x)
echo $x                  (Should be 12)

x = (4 * ($x * 10 / (5*2)) / 4) + 3
echo $x                  (Should be 15)
```
This input file can be found in the examples folder under the name arithmetic.

#### Reading file data
Data can be loaded from a file via the following syntax:
```
load <filepath> <file_type> as <variable name>

            or

read <filepath> <file_type> as <variable name>
```
Here the filepath points towards the file you would like to load, file_type is the type of file you would like to load (see example for more info) and variable name is the name associated to the data.

Many of the example input files load data, have a look at them to see more examples of the syntax.

### Setting more complex metadata
More complex metadata can be set with json files in the src/data/<set_type> folder. For example, if the molecular system you would like to analyse is pentacene you can set pentacene metadata to any data you have just loaded.

This can be achieved with the following syntax:
```
set <set_type> <variable_name> to <set_data>
```

In the pentacene example above the syntax: `set system data to pentacene` would allow the variable `data` to access pentacene metadata such as number of atoms per molecule, molecular mass, atom types etc...

To add new metadata files save them to `src/data/<set_type>/<set_data>.json`. For the pentacene example above we can save a file named `pentacene.json` in the folder `src/data/systems`. Use a json format for the data, have a look at the pentacene file to check how this should look.

An example is given in the `examples/calc_density` folder.

#### Writing file data
Data can be written to a file via the following syntax:
```
write <variable name> <filepath>

            or

write <variable name> <filepath> as <file type>
```
In the above variable name refers to the name assigned to the data to be written. Filepath is the path pointing to the file to be written and file type is the type of the file to be written.

If the first syntax is chosen the data will be written in its default file type.

Many of the example input files write data, have a look at them to see more examples of the syntax.

#### Calculating values
To calculate properties of the data the following syntax can be used.
```
calc <value> from <variable name> as <new variable name>
```
In the above line a value is being calculated from some data that has name 'variable name' and then this is saved as 'new variable name'. This new variable can overwrite the old variable.

An example of calculating something from some data is given in the calc_density folder.

#### Print to console
To print something to the console use the following syntax:
```
echo <string>
```
You can print variables in this by calling them, as in bash, with a dollar sign like:
```
echo "Value of variable: $variable"
```

#### Interacting with data with IPython
You can interact with the data you've created through the IPython shell if you use the command:
```
shell
```
Here any variables you declare in the input file (including loaded data) will be saved with a global scope and you can access them as you normally would. The shell
will open with the code's current state and after closing it will remember the altered state and run the rest of the input script with this.

An example is given in the `examples/ipython_shell` directory.
