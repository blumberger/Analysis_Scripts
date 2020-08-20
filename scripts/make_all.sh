#!/usr/bin/env bash

###############################
# This little script will make any larger modules that need making.

# It must be called from the ./run.sh script (as this export variables to it)
# It will move all the compiled executables into the bin directory

# To add new makefiles just add values to the arrays below and remember to
###############################

MISC_C_MOD_DIR="$PWD/src/misc_C_modules"
BIN_DIR="$PWD/bin"
# Declare the programs that need compiling and the necessary compilation flags

# Don't put the .c extension on this
declare -a C_DIR=("$MISC_C_MOD_DIR/AOM_Couplings")
declare -a MAKE_CMD=("all")
declare -a EXE_NAMES=("STO_proj_AOM_overlap")

NPROGS="${#C_DIR[@]}"

source $CONFIG_VAR_FILE

# Create the bin directory (if it doesn't exist)
if ! [ -d "$BIN_DIR" ]
then
    mkdir $BIN_DIR
fi

# Loop over all programs and compile 
for (( i=0; i<$NPROGS; i++ ));
do
    cd ${C_DIR[$i]};
    echo "Making ${EXE_NAMES[$i]}"
    make ${MAKE_CMD[$i]} -s
    cp $EXE_NAMES $BIN_DIR;
    cd -;
done

