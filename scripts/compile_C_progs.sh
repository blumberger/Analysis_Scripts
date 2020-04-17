#!/usr/bin/env bash

###############################
# This little script will compile any C programs used within this package.

# It must be called from the ./run.sh script (as this export variables to it)
# It will compile everything into the bin directory.
# It will only compile things if the C file has changed since the last time
#  it was checked.
# It will only compile single C files and won't create object files.


# To add new C programs just add values to the arrays below
###############################

MISC_C_MOD_DIR="$PWD/src/misc_C_modules"
BIN_DIR="$PWD/bin"
# Declare the programs that need compiling and the necessary compilation flags

# Don't put the .c extension on this
declare -a C_DIR=("$MISC_C_MOD_DIR/RDF")
declare -a C_PROGRAMS=("RDF2_inter_tri")
declare -a LIBS=("-lm")
declare -a INCLUDES=("")
declare -a CFLAGS=("-w -O3")
declare -a EXE_NAME=("RDF_inter")

NPROGS="${#C_PROGRAMS[@]}"

source $CONFIG_VAR_FILE

# Create the bin directory (if it doesn't exist)
if ! [ -d "$BIN_DIR" ]
then
    mkdir $BIN_DIR
fi

# Loop over all programs and compile 
for (( i=0; i<$NPROGS; i++ ));
do
   exe_file="$BIN_DIR/${EXE_NAME[$i]}"
   C_file="${C_DIR[$i]}/${C_PROGRAMS[$i]}.c"
   last_change_var_name="${C_PROGRAMS[$i]}"
   last_change=$(stat $C_file --printf=%Z)

   compile_cmd="${CFLAGS[$i]} ${INCLUDES[$i]} -o $exe_file $C_file ${LIBS[$i]}"

   if ! [ -f $exe_file ]
   then
      if gcc $compile_cmd;
      then
         echo "Compiled ${C_PROGRAMS[$i]}" 
         sed -i s/"$last_change_var_name=.*"/"$last_change_var_name=$last_change"/ $CONFIG_VAR_FILE
      else
         exit_code=1
         export exit_code
         exit
      fi
   else
      # Check if the last compile time is available to us
      is_there=`grep -E "$last_change_var_name=[0-9]+" $CONFIG_VAR_FILE`

      # If the last change of the C file hasn't been recorded then record it
      if [ -z "$is_there" ]
      then
         echo "$last_change_var_name=$last_change" >> $CONFIG_VAR_FILE
      fi

      if [ ! -z "$is_there" ]
      then
         last_change_prev="${!last_change_var_name}"
         if [ "$last_change_prev" != "$last_change" ]
         then
            if gcc $compile_cmd;
            then
               echo "Compiled ${C_PROGRAMS[$i]}" 
               sed -i s/"$last_change_var_name=.*"/"$last_change_var_name=$last_change"/ $CONFIG_VAR_FILE
            else
               exit_code=1
               export exit_code
               exit
            fi
         fi
      fi
   fi

done

