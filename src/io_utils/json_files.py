#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A module to handle json reading and writing
"""
import json

def write_json(data, filepath):
    """
    Will write the nearest neighbour list as a json file.
    
    Inputs:
        * data <dict> => A dict object to write to file.
        * filepath <str> => The filepath to save the data to.
    """
    # A hack to get the data from a variable
    data = dict(eval(str(data)))

    with open(filepath, 'w') as f:
        json.dump(data, f)


def read_json(filepath):
    """
    Will write the nearest neighbour list as a json file.
    
    Inputs:
        * filepath <str> => The filepath to save the data to.
    """
    with open(filepath, 'r') as f:
   	    data = json.read_json(f) 
    return data
