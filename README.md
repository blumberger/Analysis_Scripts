# MD analysis scripts

## Mission Statement
The purpose of this repo is to centralise all the scripts that we use as a group so each time someone runs an MD simulation/SH simulation they don't have to re-write what many people have written before them.

One of the main reasons the analysis scripts are re-written is (I think) because it is often easier (or less frustrating) to write the script yourself rather than learn how to use the script that's already available. Therefore, for this repo I think we need some rules about submitting that make the code both usable for people in the future and readable for anyone wanting to edit and extend the scripts.

As Guido von Russof says 'Code it read much more than it is written' so reability is important!



## The rules

 - **Use an input file**.  We are all used to using input files to interact with the various programs we use. They make things clearer for the user when it comes to running things and they can be used to standardise input across many different scripts.
 - **Use clear variable names**. Using clear variable names helps a lot when trying to understand what the script is doing. They can often replace comments (which saves time and effort for the developer)! See [PEP8](https://realpython.com/python-pep8/#naming-styles) for more info. P.S. this extends to function names, filenames, class names etc...
 - **Write with modularity in mind**. Are you writing some code that plots data from an xyz file? Then put the xyz file functions in a separate file to the plotting functions. That way people can re-use the xyz functions without having to plot a graph! It also helps to have different files/modules in sensible places! Don't put all your code in 1 file called a module called utils.py when it is very hard to know what is in utils.py. Folders can help subdivide code further too i.e. if you have 3 different modules that load xyz, pdb and .dat files these can go in a folder named 'io' (input-output).
 - **If the code can be reused, write a function**. If someone wants to write an xyz file and the function that writes them is incorporated into another function they have to do some tedious editting. Remember to always give that function a docstring though! They really help.
 - **Follow PEP8**. If you are writing in python, try and follow PEP8 guidelines as closely as possible. Most IDEs have PEP8 style checking, turn it on and try to follow it.
 - **Document your code**. It is very helpful to document your code, without the documentation people won't know if your code even exists and they will write their own! I have written a separate section on documentation below.


## How to document your code.
Some of the documentation will be automatically created through the default variables that are used in the input file. Though a lot of it can't be created with a nice script i.e. how the code works, how to edit it and what it is used for. For this one should create a folder within the directory 'examples' with markdown files contained inside. These markdown files should explain how the code works, how to edit it and give and example of how it would be ran.

## What should go in the repo.
Any code you think may be helpful, can be ran with an input file and is documented.

