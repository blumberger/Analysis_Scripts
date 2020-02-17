CONFIG_DIR="./config"
CONFIG_VAR_FILE="$CONFIG_DIR/vars"
PWD=`pwd`


# Try to read the variable config file, this is just a file that stores variables permanently. It will be ignored by git.
INSTALL_DEPS="false"
if ! [ -f "$CONFIG_VAR_FILE" ]
then
    INSTALL_DEPS="true"
else
    source $CONFIG_VAR_FILE
fi

# Override the config file if pipenv isn't installed
PIPENV_EXISTS=`which pipenv 2> /dev/null`
if [ "$PIPENV_EXISTS" == "" ]
then
   INSTALL_DEPS="true"
fi

# Install the necessary python libraries and create a virtual enviroment
if [ "$INSTALL_DEPS" == "true" ]
then
    echo "Installing dependencies"
    if [ "$PIPENV_EXISTS" == "" ]
    then
        python3 -m pip install pipenv --user
    fi
    pipenv install
    echo "Installed dependencies"
fi

# Create the config directory if it doesn't exist
if ! [ -d "$CONFIG_DIR" ]
then
   mkdir $CONFIG_DIR
fi

# Let the script know we don't need to install things next time.
echo "INSTALL_DEPS=\"false\"" &> $CONFIG_VAR_FILE
    
# Get the input file from the arguments
INP_FILE=""
TEST="false"
HELP="true"
while getopts "i:thu" arg
do
    case $arg in 
        i) INP_FILE=$OPTARG; HELP="false";;
        t) TEST="true"; HELP="false";;
        u) HELP="false"; pipenv install;;
        h) HELP="true";;
    esac
done

if [ "$HELP" == "true" ]
then
    echo " "
    echo " "
    echo "/---------- MD Analysis Utility Help ------------------\\"
    echo "|                                                      |"
    echo "| To use this utility please supply an input file.     |"
    echo "|                                                      |"
    echo "| To find get some inspiration or to check the syntax  |"
    echo "| of input files check out the examples folder. A full |"
    echo "| list of input foldernames is below:                  |"
    echo "|                                                      |"
    python3 -c "s=\"\"\"$(ls examples/)\"\"\"; l=[f'{j})  {i}'.ljust(47)+'|' for j, i in enumerate(s.split()) if 'data' not in i]; print(\"|\t\" + '\n|\t'.join(l))"
    echo "|                                                      |"
    echo "| These input files have been well commented and       |"
    echo "| all operations are explained.                        |"
    echo "|                                                      |"
    echo "|------------------------------------------------------|"
    echo "|                                                      |"
    echo "| Running the Code.                                    |"
    echo "|------------------------------------------------------|"
    echo "| In order to run the code simply run this run.sh file |"
    echo "| with bash, specifying where the input file is with   |"
    echo "| the flag '-i'.                                       |"
    echo "|                                                      |"
    echo "| E.g:     './run.sh -i <input_filepath>'              |"
    echo "|                                                      |"
    echo "|                                                      |"
    echo "|------------------------------------------------------|"
    echo "|                                                      |"
    echo "| Flags                                                |"
    echo "|------------------------------------------------------|"
    echo "| -t  -> Test the code                                 |"
    echo "| -i  -> Specify input file filepath                   |"
    echo "| -u  -> Update dependencies                           |"
    echo "|                                                      |"
    echo "|------------------------------------------------------|"
    echo "|                                                      |"
    echo "| Getting more help.                                   |"
    echo "|------------------------------------------------------|"
    echo "| If you would like more info read the README          |"
    echo "| or ask Matt                                          |"
    echo "\\------------------------------------------------------/"
    exit 1
fi


if [ "$TEST" == "true" ]
then
    pipenv run python3 $PWD/src/tests/call_all_examples.py
    exit 0
fi


# If we can't find the input file raise an error
if [ "$INP_FILE" == "" ]
then
    exit 0
fi
if ! [ -f "$INP_FILE" ]
then
    echo "ERROR: Can't find file '"$INP_FILE"'"
    exit 2
fi

# If everything is ok, pass the input file to the python code
pipenv run python3 $PWD/main.py -i $INP_FILE
