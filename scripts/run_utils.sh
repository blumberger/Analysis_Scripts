# Will install any dependencies and init the pip virtual environment
function install_deps() {
    echo "Installing dependencies"
    if [ "$PIPENV_EXISTS" == "" ]
    then
        python3 -m pip install pipenv --user
    fi  
    pipenv install
    pipenv install ipython
    echo "Installed dependencies"
}

