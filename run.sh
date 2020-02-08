CONFIG_VAR_FILE="./config/vars"


# Try to read the variable config file, this is just a file that stores variables permanently. It will be ignored by git.
INSTALL_DEPS="false"
if ! [ -f $CONFIG_VAR_FILE ]
then
    INSTALL_DEPS="true"
else
    source $CONFIG_VAR_FILE
fi

# Install the necessary python libraries and create a virtual enviroment
if [ $INSTALL_DEPS == "true" ]
then
    echo "Installing dependencies"
    PIPENV_EXISTS=`which pipenv`
    if [ $PIPENV_EXISTS == "" ]
    then
        python3 -m pip install pipenv --user
    fi
    pipenv install
    echo "Installed dependencies"
fi
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
