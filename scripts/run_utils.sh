# Will install any python dependencies and init the pip virtual environment
function install_deps() {
    echo "Installing dependencies"
    
    exit_code=0

    if [ "$PIPENV_EXISTS" == "" ]
    then
        python3 -m pip install pipenv --user --force-reinstall --no-cache-dir --upgrade
    fi  
    
    # If the pipenv exe still doesn't exist
    PIPENV_EXISTS=`which pipenv 2> /dev/null`
    if [ "$PIPENV_EXISTS" == "" ]
    then
        # Check the place the executable is likely to be in
        if [ -f "$HOME/.local/bin/pipenv" ]
        then
            # Add the path to the bashrc file and source it
            echo "Found your pipenv executable in '$HOME/.local/bin/pipenv'"
            echo "Adding this as part of your path in the ~/.bashrc file"

            echo "######################################################" >> ./bashrc
            echo " " >> ./bashrc
            echo "######################################################" >> ./bashrc
            echo "# Adding for python binaries (pipenv and pip etc...) #" >> ./bashrc
            echo "######################################################" >> ./bashrc
            echo "PATH=\"$HOME/.local/bin/pipenv:\$PATH\"" >> ~/.bashrc
            echo "export PATH" >> ~/.bashrc
            echo "######################################################" >> ./bashrc
            echo " " >> ./bashrc

        else
            # Raise an error
            echo "Can't find the pipenv executable."
            echo "It doesn't seem to have been installed through the command 'python3 -m pip install pipenv --user --force-reinstall --no-cache-dir --upgrade'"
            echo "Try installing it yourself to use the script"
            exit_code=1
        fi
    fi
    
    pipenv install ipython
    pipenv install
    echo "Installed dependencies"

    export exit_code 
}
