#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains all the code relevant to setting the system parameters in
Variable objects.
"""
import os

from src.io_utils import json_files

SYSTEMS_FOLDER = "src/data/systems"
SYSTEMS_FILES = os.listdir(SYSTEMS_FOLDER)

def set_system(Var, sys_name):
    """
    Will set system data to the metadata in a Variable.

    Inputs:
        * Var <Variable> => An instance of the Variable class
        * sys_name <str> => The string containing the system to get data from.
    Outputs:
        <Variable>, <int> The variable with changed properties and exit code
    """
    sys_name = sys_name.replace(" ", "_").lower() + ".json"

    # Filepath doesn't exist
    if sys_name not in SYSTEMS_FILES:
        return False, 1

    system_filepath = os.path.join(SYSTEMS_FOLDER, sys_name)
    system = json_files.read_json(system_filepath)
    for key in system:
        Var.metadata[key] = system[key]

    return Var, 0
