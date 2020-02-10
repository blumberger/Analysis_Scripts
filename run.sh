CONFIG_DIR="./config"
CONFIG_VAR_FILE="$CONFIG_DIR/vars"


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
while getopts "i:" arg
do
    case $arg in 
        i) INP_FILE=$OPTARG;;
    esac
done

# If we can't find the input file raise an error
if [ "$INP_FILE" == "" ]
then
    echo "Please supply an input file to the run.sh as \`./run.sh -i "input.inp"\`"
    exit 1
fi
if ! [ -f $INP_FILE ]
then
    echo "ERROR: Can't find file '"$INP_FILE"'"
    exit 2
fi

# If everything is ok, pass the input file to the python code
pipenv run python3 main.py -i $INP_FILE
